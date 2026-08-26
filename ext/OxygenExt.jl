# ext/OxygenExt.jl -- kora-server: Oxygen.jl HTTP application around Kora.jl.
#
# Weakdep extension (matches ArrowExt/MakieExt) so Oxygen/HTTP only enter the
# dependency tree of consumers that `using Oxygen`; the minimal kora-worker
# binary and the native desktop bridge are unaffected.
#
# Routes, per .claude/plans/web-app/kora-web-service.md Component 1:
#   POST /api/session/start -> { session_token }
#   POST /api/run_reef      -> WireSimParams bytes in, WireEnsembleResult bytes out
#   POST /api/session/end   -> { ok: true }
#   GET  /health             -> { status, version }
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

# ---------------------------------------------------------------------------
# Wire sizes -- must stay in sync with build/worker_main.jl
# ---------------------------------------------------------------------------
const WORKER_PARAMS_BYTES = 48
const WORKER_RESULT_BYTES = 34804

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
    return HTTP.Response(
        status, ["Content-Type" => "application/json; charset=utf-8"],
        body="{\"error\": \"$msg\"}"
    )
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

    cors = Oxygen.Cors(; allowed_origins=cors_allowed_origins)

    Oxygen.serve(; host=host, port=port, middleware=[cors], async=false)
end

end  # module OxygenExt
