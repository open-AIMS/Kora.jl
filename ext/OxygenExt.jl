# ext/OxygenExt.jl -- kora-server: Oxygen.jl HTTP application around Kora.jl.
#
# Weakdep extension (matches ArrowExt/MakieExt) so Oxygen/HTTP only enter the
# dependency tree of consumers that `using Oxygen`; the minimal kora-worker
# binary and the native desktop bridge are unaffected.
#
# Routes, per .claude/plans/web-app/kora-web-service.md Component 1:
#   POST /api/session/start   -> { session_token }
#   POST /api/run_reef        -> WireSimParams bytes in, WireEnsembleResult bytes out
#   POST /api/session/end     -> { ok: true }
#   POST /api/dhw/parse_netcdf -> raw NetCDF/HDF5 bytes in, parsed DHW cube as JSON out
#                                 (Part 3's wasm-NetCDF increment -- see below)
#   GET  /health               -> { status, version }
#
# kora-server spawns one kora-worker subprocess per session (heap isolation;
# a crash in one user's simulation cannot affect others) and talks to it over
# the child's stdin/stdout using the fixed-size binary framing worker_main.jl
# already implements -- no length prefix needed since both message shapes are
# constant size.
module OxygenExt

using Kora
using Oxygen
using HTTP
using UUIDs
using NCDatasets

# ---------------------------------------------------------------------------
# Wire sizes -- must stay in sync with build/worker_main.jl
# ---------------------------------------------------------------------------
const WORKER_PARAMS_BYTES = 1420  # Part 3: +dhw_override[75] f32 +dhw_override_active u32
                                   # Part 5 v2: +init_size_class_fraction[35] f32 +init_size_class_active u32
                                   # Part 10: dhw_override grown 75->MAX_TIMESTEPS(300) f32, +n_timesteps u32
const WORKER_RESULT_BYTES = 157208  # Part 10: every time-indexed field grown to MAX_TIMESTEPS(300)
                                     # capacity, +trailing n_timesteps u32 -- see worker_main.jl's header

# ---------------------------------------------------------------------------
# Session / worker-process registry
# ---------------------------------------------------------------------------
mutable struct WorkerHandle
    proc::Base.Process
    io_lock::ReentrantLock
    last_used::Float64
end

const _WORKERS = Dict{String,WorkerHandle}()
const _WORKERS_LOCK = ReentrantLock()

# ---------------------------------------------------------------------------
# Config resolution -- environment variables, with sensible local-dev
# fallbacks. ReefGuide-specific defaults (production CORS origin, Sentry,
# etc.) are KoraReefGuideWorker.jl's responsibility, not kora-server's --
# see Component 1 in the plan.
# ---------------------------------------------------------------------------
function _default_worker_bin()::String
    env = get(ENV, "KORA_WORKER_BIN", "")
    isempty(env) || return env

    # juliac --bundle deployment layout: kora-server and kora-worker are
    # copied side-by-side into the same bin/ dir (see Phase 4 Dockerfile).
    try
        self_exe = readlink("/proc/self/exe")
        sibling = joinpath(dirname(self_exe), "kora-worker")
        isfile(sibling) && return sibling
    catch
        # not on Linux, or /proc unavailable (e.g. plain `julia` REPL) -- fall through
    end

    # Local dev fallback: the worker binary produced by `build.sh --mode worker`.
    dev_path = joinpath(pkgdir(Kora), "build", "dist", "worker", "bin", "kora-worker")
    isfile(dev_path) && return dev_path

    return "kora-worker"  # last resort: rely on PATH
end

function _default_model_path(envvar::String, filename::String)::String
    env = get(ENV, envvar, "")
    isempty(env) || return env
    return joinpath(Kora._kora_assets_dir(), "models", filename)
end

_default_growth_model_path() =
    _default_model_path("KORA_GROWTH_MODEL_PATH", "offshore_north_growth_models.json")
_default_survival_model_path() =
    _default_model_path("KORA_SURVIVAL_MODEL_PATH", "offshore_north_survival_models.json")

_env_int(name::String, default::Int)::Int = parse(Int, get(ENV, name, string(default)))

# ---------------------------------------------------------------------------
# Worker process lifecycle
# ---------------------------------------------------------------------------
function _spawn_worker(worker_bin::String, growth_path::String, survival_path::String)::WorkerHandle
    proc = open(`$worker_bin $growth_path $survival_path`, "r+")
    return WorkerHandle(proc, ReentrantLock(), time())
end

function _kill_worker!(handle::WorkerHandle)::Nothing
    try
        close(handle.proc.in)
    catch
    end
    try
        kill(handle.proc)
    catch
    end
    return nothing
end

# Read exactly n bytes, blocking as needed; throws EOFError if the worker
# exits early (e.g. crashed mid-simulation).
function _read_exact!(io::IO, n::Int)::Vector{UInt8}
    buf = Vector{UInt8}(undef, n)
    read!(io, buf)
    return buf
end

function _run_on_worker(handle::WorkerHandle, params_bytes::Vector{UInt8})::Vector{UInt8}
    length(params_bytes) == WORKER_PARAMS_BYTES || error(
        "params size mismatch: got $(length(params_bytes)), expected $WORKER_PARAMS_BYTES"
    )
    lock(handle.io_lock) do
        write(handle.proc, params_bytes)
        flush(handle.proc)
        _read_exact!(handle.proc, WORKER_RESULT_BYTES)
    end
end

# ---------------------------------------------------------------------------
# Idle sweep -- kills worker processes unused for IDLE_TIMEOUT_MS (default
# 2 min). This is the per-worker timeout described in the plan's
# "kora-server lifecycle & cost" section; kora-server's OWN idle/scale
# lifecycle is a separate, KoraReefGuideWorker.jl-level concern (Phase 4/6).
# ---------------------------------------------------------------------------
function _sweep_idle_workers!(idle_timeout_s::Float64)::Nothing
    now = time()
    lock(_WORKERS_LOCK) do
        for (token, handle) in collect(_WORKERS)
            if now - handle.last_used > idle_timeout_s
                _kill_worker!(handle)
                delete!(_WORKERS, token)
            end
        end
    end
    return nothing
end

function _start_idle_sweeper(idle_timeout_s::Float64)::Timer
    # Sweep at 1/4 the idle timeout, floored at 5s, so the actual time a dead
    # session lingers is bounded close to idle_timeout_s.
    interval = max(5.0, idle_timeout_s / 4)
    return Timer(interval; interval=interval) do _
        _sweep_idle_workers!(idle_timeout_s)
    end
end

# ---------------------------------------------------------------------------
# Auth helper -- session token carried as `Authorization: Bearer <token>`
# ---------------------------------------------------------------------------
function _bearer_token(req::HTTP.Request)::Union{String,Nothing}
    auth = HTTP.header(req, "Authorization", "")
    startswith(auth, "Bearer ") || return nothing
    token = auth[8:end]
    isempty(token) && return nothing
    return token
end

function _json_error(status::Int, msg::String)::HTTP.Response
    # Every prior call site passed a fixed short literal (no quotes/control
    # chars), so plain interpolation never mattered -- but /api/dhw/parse_netcdf
    # forwards a `showerror` message that can contain anything (a quoted file
    # path, an embedded newline), so escape properly here.
    escaped = replace(msg, "\\" => "\\\\", "\"" => "\\\"", "\n" => "\\n", "\r" => "\\r")
    return HTTP.Response(
        status, ["Content-Type" => "application/json; charset=utf-8"],
        body="{\"error\": \"$escaped\"}"
    )
end

# ---------------------------------------------------------------------------
# /api/dhw/parse_netcdf -- Part 3's wasm-NetCDF increment. wasm can't link
# libnetcdf directly (see kora-app's `dhw_netcdf.rs` header comment), so this
# stateless endpoint does the parsing server-side with NCDatasets (the server
# exe is untrimmed, so this is viable here -- see the plan doc's OQ-2
# resolution) and hands back the (time, location) matrix as JSON. No session
# token needed: this doesn't touch a kora-worker process.
#
# Dimension classification mirrors `dhw_netcdf.rs` exactly -- by *name*, never
# position, since real ADRIA cubes declare `dhw(member, location, timesteps)`
# on disk (see that file's header comment for why). NCDatasets already
# presents `dimnames`/indexing in Julia's own (column-major) convention, so
# unlike the Rust reader (which reads the raw C-order buffer by hand) no
# manual byte-order reversal is needed here -- just axis bookkeeping.
# ---------------------------------------------------------------------------
const _DHW_TIMESTEP_DIM_NAMES = ("time", "timestep", "timesteps", "year", "years", "t")
const _DHW_MEMBER_DIM_PATTERNS = ("member", "scenario", "draw", "ensemble", "realisation", "realization")
const _DHW_LOCATION_LABEL_PREFERRED =
    ("reef_siteid", "site_id", "siteid", "location", "locations", "site", "unique_id", "name", "id")

_is_dhw_timestep_dim(name::AbstractString) = lowercase(name) in _DHW_TIMESTEP_DIM_NAMES
_is_dhw_member_dim(name::AbstractString) =
    any(p -> occursin(p, lowercase(name)), _DHW_MEMBER_DIM_PATTERNS)

# The `dhw` (or `DHW`) variable, or -- failing that -- the sole floating-point
# variable with 2 or 3 dimensions. Mirrors `find_dhw_variable` in dhw_netcdf.rs.
function _find_dhw_variable(ds::NCDataset)::String
    for name in ("dhw", "DHW")
        haskey(ds, name) && return name
    end
    candidates = String[]
    for name in keys(ds)
        v = ds[name]
        et = eltype(v)
        et2 = et isa Union ? Base.uniontypes(et) : (et,)
        is_float = any(t -> t <: AbstractFloat, et2)
        if is_float && ndims(v) in (2, 3)
            push!(candidates, name)
        end
    end
    if length(candidates) == 1
        return candidates[1]
    elseif isempty(candidates)
        error("no 'dhw' variable found, and no 2-D/3-D floating-point variable to fall back to")
    else
        error(
            "no 'dhw' variable found, and $(length(candidates)) candidate floating-point " *
            "variables are ambiguous (expected exactly one)"
        )
    end
end

# Location labels from whichever string-valued variable shares exactly the
# location dimension, preferring the ADRIA-conventional names. `nothing` when
# no such variable exists, or a read fails partway (a partial label set would
# be more confusing than none) -- mirrors `find_location_labels`.
function _find_dhw_location_labels(
    ds::NCDataset, loc_dim_name::String, n_loc::Int
)::Union{Vector{String},Nothing}
    candidates = String[]
    for name in keys(ds)
        v = ds[name]
        if dimnames(v) == (loc_dim_name,) && eltype(v) <: Union{AbstractString,Missing}
            push!(candidates, name)
        end
    end
    isempty(candidates) && return nothing
    sort!(candidates; by=name -> something(
        findfirst(==(lowercase(name)), _DHW_LOCATION_LABEL_PREFERRED),
        length(_DHW_LOCATION_LABEL_PREFERRED) + 1,
    ))
    try
        raw = ds[candidates[1]][:]
        length(raw) == n_loc || return nothing
        return [ismissing(x) ? "" : String(x) for x in raw]
    catch
        return nothing
    end
end

# Parse a NetCDF/HDF5 DHW trajectory from raw bytes. Returns a `Dict` matching
# `NetcdfDhwParseResponse` on the Rust side: `n_timesteps`, `n_locations`,
# `location_labels` (`Vector{String}` or `nothing`), `values` (row-major
# `(time, location)`, i.e. `values[t * n_locations + loc]`, 0-based, exactly
# like `dhw_csv::parse_csv`/`dhw_netcdf::parse_netcdf`). Throws on any parse
# error -- the caller turns that into a 400.
function _parse_dhw_netcdf_bytes(bytes::Vector{UInt8})::Dict{String,Any}
    tmppath = tempname() * ".nc"
    write(tmppath, bytes)
    try
        NCDataset(tmppath) do ds
            varname = _find_dhw_variable(ds)
            var = ds[varname]
            dnames = collect(dimnames(var))
            nd = length(dnames)
            if nd != 2 && nd != 3
                error(
                    "'$varname' has $nd dimensions ($(join(dnames, ", "))); expected 2 " *
                    "(time, location) or 3 (time, location, member)"
                )
            end

            time_idx = findfirst(_is_dhw_timestep_dim, dnames)
            time_idx === nothing && error(
                "could not find a timestep dimension on '$varname' (dims: $(join(dnames, ", ")))"
            )

            other = [i for i in 1:nd if i != time_idx]
            local loc_idx, member_idx
            if length(other) == 1
                loc_idx, member_idx = other[1], nothing
            else
                mpos = findfirst(i -> _is_dhw_member_dim(dnames[i]), other)
                mpos === nothing && error(
                    "'$varname' has 3 dimensions but none of them look like a " *
                    "member/scenario axis (expected a name containing 'member', " *
                    "'scenario', 'draw' or 'ensemble'; dims: $(join(dnames, ", ")))"
                )
                member_idx = other[mpos]
                loc_idx = other[mpos == 1 ? 2 : 1]
            end

            n_time = size(var, time_idx)
            n_loc = size(var, loc_idx)
            loc_dim_name = dnames[loc_idx]

            idx = Any[Colon() for _ in 1:nd]
            member_idx === nothing || (idx[member_idx] = 1)  # select member 0 (1-based here)
            raw = Array(var[idx...])  # drops the member axis, if any

            # Dropping (at most) the member axis preserves the relative order
            # of the remaining two axes, so the pre-drop time_idx/loc_idx
            # comparison still tells us whether `raw` is (time, location) or
            # (location, time).
            canon = Float32.(time_idx < loc_idx ? raw : permutedims(raw, (2, 1)))  # now (time, location)
            values = [canon[t, loc] for t in 1:n_time for loc in 1:n_loc]  # row-major flatten

            labels = _find_dhw_location_labels(ds, loc_dim_name, n_loc)

            return Dict{String,Any}(
                "n_timesteps" => n_time,
                "n_locations" => n_loc,
                "location_labels" => labels,
                "values" => values,
            )
        end
    finally
        rm(tmppath; force=true)
    end
end

# ---------------------------------------------------------------------------
# Kora.start_server -- the public entry point, called from a plain `julia`
# session (`using Kora, Oxygen; Kora.start_server()`) or from the
# juliac-compiled build/kora_server_main.jl entry point.
# ---------------------------------------------------------------------------
function Kora.start_server(;
    host::String=get(ENV, "KORA_SERVER_HOST", "0.0.0.0"),
    port::Int=_env_int("KORA_SERVER_PORT", 4444),
    worker_bin::String=_default_worker_bin(),
    growth_model_path::String=_default_growth_model_path(),
    survival_model_path::String=_default_survival_model_path(),
    idle_timeout_ms::Int=_env_int("IDLE_TIMEOUT_MS", 120_000),
    max_concurrent_workers::Int=_env_int("MAX_CONCURRENT_WORKERS", 20),
    cors_allowed_origins::Vector{String}=String.(split(get(ENV, "CORS_ALLOWED_ORIGINS", "*"), ",")),
)
    idle_timeout_s = idle_timeout_ms / 1000

    @get "/health" function(req::HTTP.Request)
        return Dict("status" => "ok", "version" => string(pkgversion(Kora)))
    end

    @post "/api/session/start" function(req::HTTP.Request)
        n_active = lock(() -> length(_WORKERS), _WORKERS_LOCK)
        if n_active >= max_concurrent_workers
            return _json_error(503, "server_at_capacity")
        end

        handle = try
            _spawn_worker(worker_bin, growth_model_path, survival_model_path)
        catch e
            @error "failed to spawn kora-worker" exception = (e, catch_backtrace())
            return _json_error(500, "worker_spawn_failed")
        end

        token = string(uuid4())
        lock(_WORKERS_LOCK) do
            _WORKERS[token] = handle
        end
        return Dict("session_token" => token)
    end

    @post "/api/run_reef" function(req::HTTP.Request)
        token = _bearer_token(req)
        token === nothing && return _json_error(401, "missing_session_token")

        handle = lock(() -> get(_WORKERS, token, nothing), _WORKERS_LOCK)
        handle === nothing && return _json_error(404, "unknown_session")

        params_bytes = Oxygen.binary(req)
        if params_bytes === nothing || length(params_bytes) != WORKER_PARAMS_BYTES
            return _json_error(400, "bad_params_length")
        end

        result_bytes = try
            _run_on_worker(handle, params_bytes)
        catch e
            @error "kora-worker call failed" exception = (e, catch_backtrace())
            lock(_WORKERS_LOCK) do
                delete!(_WORKERS, token)
            end
            _kill_worker!(handle)
            return _json_error(502, "worker_unavailable")
        end

        handle.last_used = time()
        return HTTP.Response(
            200, ["Content-Type" => "application/octet-stream"], body=result_bytes
        )
    end

    @post "/api/dhw/parse_netcdf" function(req::HTTP.Request)
        # `Oxygen.binary` is declared to return `Vector{UInt8}` but actually
        # returns `nothing` on a genuinely empty body -- a `TypeError` on the
        # way out, not a value our `=== nothing` check below ever gets to see.
        # Same latent bug as `/api/run_reef`'s equivalent check; caught here
        # explicitly so an empty request gets a clean 400 instead of a 500.
        bytes = try
            Oxygen.binary(req)
        catch
            nothing
        end
        if bytes === nothing || isempty(bytes)
            return _json_error(400, "empty_body")
        end
        result = try
            _parse_dhw_netcdf_bytes(bytes)
        catch e
            msg = sprint(showerror, e)
            return _json_error(400, "parse_failed: $msg")
        end
        return result
    end

    @post "/api/session/end" function(req::HTTP.Request)
        token = _bearer_token(req)
        token === nothing && return _json_error(401, "missing_session_token")

        handle = lock(_WORKERS_LOCK) do
            pop!(_WORKERS, token, nothing)
        end
        handle !== nothing && _kill_worker!(handle)
        return Dict("ok" => true)
    end

    _start_idle_sweeper(idle_timeout_s)

    # Oxygen.Cors' `allowed_headers` defaults to ["*"], but the CORS spec
    # carves out an exception for wildcards: "*" never authorizes the
    # `Authorization` header (fetch spec 3.2.5 / whatwg). run_reef and
    # session/end send `Authorization: Bearer <token>`, so browsers silently
    # block those requests at the preflight stage with a bare CORS failure
    # (Firefox: "NetworkError when attempting to fetch resource"; Chrome:
    # "Failed to fetch") even though curl/non-browser clients work fine.
    # List the actual headers instead of relying on the wildcard.
    cors = Oxygen.Cors(;
        allowed_origins=cors_allowed_origins,
        allowed_headers=["Authorization", "Content-Type"],
    )

    Oxygen.serve(; host=host, port=port, middleware=[cors], async=false)
end

end  # module OxygenExt
