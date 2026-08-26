# worker_main.jl -- juliac --output-exe entry point for kora-worker
#
# Protocol (stdin/stdout, binary):
#   For each simulation request:
#     1. Coordinator writes WORKER_PARAMS_BYTES bytes of WorkerSimParams to stdin.
#     2. Worker reads exactly WORKER_PARAMS_BYTES bytes, runs the simulation,
#        writes exactly WORKER_RESULT_BYTES bytes of WorkerEnsembleResult to stdout.
#   Loop until stdin closes.
#
# All text diagnostics go to stderr; stdout is purely binary.
#
# WorkerSimParams wire layout (little-endian, 48 bytes, no padding):
#   reef_area_m2:         f32  offset  0
#   init_cover_pct:       f32  offset  4
#   deploy_volumes[5]:    u32  offset  8  (20 bytes)
#   deploy_start_year:    u32  offset 28
#   deploy_cadence_years: u32  offset 32
#   depth_m:              u32  offset 36
#   deploy_dhw_tolerance: f32  offset 40
#   dhw_seed:             u32  offset 44
#   Total: 48 bytes
#
# NOTE: WorkerSimParams is a superset of sim-types/src/wire.rs WireSimParams —
# it adds depth_m and deploy_dhw_tolerance. wire.rs's WireSimParams matches
# this layout as of kora-app's feat/web-backend branch.
#
# WorkerEnsembleResult wire layout (little-endian, 34804 bytes, no padding):
#   n_valid_runs:                     u32  offset      0  (4 bytes)
#   covers[MAX_RUNS * N_TIMESTEPS]:   f32  offset      4  (30000 bytes), run-major
#   summary.lower[N_TIMESTEPS][N_GROUPS]:  f32  offset  30004  (1500 bytes)
#   summary.median[N_TIMESTEPS][N_GROUPS]: f32  offset  31504  (1500 bytes)
#   summary.upper[N_TIMESTEPS][N_GROUPS]:  f32  offset  33004  (1500 bytes)
#   dhw[N_TIMESTEPS]:                 f32  offset  34504  (300 bytes)
#   Total: 34804 bytes
#
# summary layout mirrors WireGroupSummary in wire.rs:
#   [[f32; N_GROUPS]; N_TIMESTEPS] = row-major with timestep as outer index.
#
# Usage:
#   kora-worker <growth_model_path> <survival_model_path>

module KoraWorker

using Kora
using Random: Xoshiro
using Statistics: quantile

# ---------------------------------------------------------------------------
# Wire layout constants — must stay in sync with sim-types/src/wire.rs
# ---------------------------------------------------------------------------
const N_GROUPS = 5
const N_TIMESTEPS = 75
const MAX_RUNS = 100

# 4 scalar fields (2x f32 + 5x u32 deploy_volumes + u32 start + u32 cadence) plus depth_m (u32), dhw_tol (f32), dhw_seed (u32)
const WORKER_PARAMS_BYTES = 4 + 4 + N_GROUPS * 4 + 4 + 4 + 4 + 4 + 4   # = 48

# u32 n_valid + [MAX_RUNS * N_TIMESTEPS] f32 covers + [N_TIMESTEPS * N_GROUPS * 3] f32 summary + [N_TIMESTEPS] f32 dhw
const WORKER_RESULT_BYTES = 4 + MAX_RUNS * N_TIMESTEPS * 4 + N_GROUPS * N_TIMESTEPS * 3 * 4 + N_TIMESTEPS * 4  # = 34804

# ---------------------------------------------------------------------------
# Global simulation state (same pattern as bridge_aot.jl)
# ---------------------------------------------------------------------------
const _growth_ref = Ref{Union{Nothing,Kora.PolyGrowthModel{Float32}}}(nothing)
const _survival_ref = Ref{Union{Nothing,Kora.PolySurvivalModel{Float32}}}(nothing)
const _dhw_ref = Ref{Union{Nothing,Matrix{Float32}}}(nothing)
const _init_n_ts_ref = Ref{Int}(0)
const _dhw_seed_ref = Ref{UInt32}(0)

# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

# libuv (which spawns this process when run as a coordinator's child, e.g.
# from kora-server) sets pipe stdio to O_NONBLOCK. A blocking-style raw
# ccall(:read) on fd 0 would then intermittently see EAGAIN (errno -11)
# before the coordinator's next write arrives -- indistinguishable from EOF
# (ret <= 0) unless we either check errno or just clear O_NONBLOCK once up
# front. The latter is simpler and keeps read_exact_stdin's "ret <= 0 means
# EOF" logic correct for both invocation styles (`< file` redirection, where
# the fd was already blocking, and a live coordinator pipe).
function _ensure_blocking_stdin()::Nothing
    F_GETFL = Cint(3)
    F_SETFL = Cint(4)
    O_NONBLOCK = Cint(0o4000)
    flags = ccall(:fcntl, Cint, (Cint, Cint), Cint(0), F_GETFL)
    flags >= 0 || return nothing
    if (flags & O_NONBLOCK) != 0
        ccall(:fcntl, Cint, (Cint, Cint, Cint), Cint(0), F_SETFL, flags & ~O_NONBLOCK)
    end
    return nothing
end

# trim-safe variant: Core.stdin resolves to Any at compile-time so Julia IO
# dispatch on it fails --trim=safe.  Use ccall(:read) on fd 0 directly instead.
function read_exact_stdin(n::Int)::Union{Vector{UInt8},Nothing}
    buf = Vector{UInt8}(undef, n)
    total = 0
    while total < n
        ret = GC.@preserve buf ccall(
            :read, Cssize_t,
            (Cint, Ptr{UInt8}, Csize_t),
            Cint(0), pointer(buf, total + 1), Csize_t(n - total)
        )
        ret <= Cssize_t(0) && return nothing  # EOF (0) or error (< 0)
        total += Int(ret)
    end
    return buf
end

# ---------------------------------------------------------------------------
# Parse WorkerSimParams from 44 raw bytes (little-endian field order)
# ---------------------------------------------------------------------------
function parse_params(bytes::Vector{UInt8})
    length(bytes) == WORKER_PARAMS_BYTES || error(
        "params size mismatch: got $(length(bytes)), expected $WORKER_PARAMS_BYTES"
    )
    io = IOBuffer(bytes)
    reef_area_m2 = read(io, Float32)
    init_cover_pct = read(io, Float32)
    deploy_volumes = ntuple(_ -> read(io, UInt32), N_GROUPS)
    deploy_start_year = read(io, UInt32)
    deploy_cadence_years = read(io, UInt32)
    depth_m = read(io, UInt32)
    deploy_dhw_tolerance = read(io, Float32)
    dhw_seed = read(io, UInt32)
    return (;
        reef_area_m2,
        init_cover_pct,
        deploy_volumes,
        deploy_start_year,
        deploy_cadence_years,
        depth_m,
        deploy_dhw_tolerance,
        dhw_seed
    )
end

# ---------------------------------------------------------------------------
# Simulation helpers (adapted from bridge_aot.jl)
# ---------------------------------------------------------------------------
function _build_ensemble_params(
    area_m2::Float32, init_cover_pct::Float32, n_members::Int
)::Matrix{Float64}
    mean_cov = Float64(Kora.mean_colony_cover_m2())
    target_cover_m2 = (Float64(init_cover_pct) / 100.0) * Float64(area_m2)
    target_pop = max(5, ceil(Int64, target_cover_m2 / mean_cov))
    pop_density = Float64(target_pop) / Float64(area_m2)
    params = Matrix{Float64}(undef, 6, n_members)
    params[1, :] .= pop_density
    params[2:6, :] .= 0.2
    return params
end

function _apply_deployment!(
    reef::Kora.ReefState, n_ts::Int,
    vols::NTuple{5,UInt32}, start::Int, cadence::Int
)::Nothing
    if start >= 1 && cadence >= 1
        for ts in start:cadence:n_ts
            for grp in 1:5
                reef.deployment_times[ts, 1, grp] = Float32(vols[grp])
            end
        end
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Simulation entry — returns WORKER_RESULT_BYTES raw bytes
# ---------------------------------------------------------------------------
function run_simulation(p)::Vector{UInt8}
    n_ts = N_TIMESTEPS

    gm = _growth_ref[]
    sm = _survival_ref[]
    (gm === nothing || sm === nothing) && error("models not loaded")

    # Reuse cached DHW unless n_ts or the requested seed changed (same policy
    # as bridge_aot.jl, plus seed-awareness -- kora-server has no separate
    # regenerate-DHW endpoint, so a client asking for a different seed on an
    # otherwise ordinary /api/run_reef call is how "New DHW trajectory" is
    # expressed here).
    if _dhw_ref[] === nothing || _init_n_ts_ref[] != n_ts || _dhw_seed_ref[] != p.dhw_seed
        _init_n_ts_ref[] = n_ts
        _dhw_seed_ref[] = p.dhw_seed
        _dhw_ref[] = Kora.generate_example_dhw(n_ts, 1; rng=Xoshiro(Int(p.dhw_seed)))
    end
    dhw_mat = _dhw_ref[]::Matrix{Float32}

    reef = Kora.initialize_reef(;
        n_timesteps=n_ts,
        n_locs=1,
        area=Float64(p.reef_area_m2),
        density=10,
        depths=Float64(p.depth_m),
        growth_models=gm,
        survival_models=sm
    )
    Kora.initialize_coral_population!(reef)
    _apply_deployment!(
        reef, n_ts, p.deploy_volumes, Int(p.deploy_start_year), Int(p.deploy_cadence_years)
    )

    n_members = 25
    ensemble_params = _build_ensemble_params(p.reef_area_m2, p.init_cover_pct, n_members)
    results = Kora.run_ensemble!(
        reef, dhw_mat, ensemble_params; deploy_dhw_tol=p.deploy_dhw_tolerance
    )

    valid_mask = [!any(isnan, results.cover[:, 1, r]) for r in 1:n_members]
    valid_indices = findall(valid_mask)
    n_valid = length(valid_indices)

    # Pre-compute per-group percentiles in one pass to avoid triple iteration.
    # Stored as [timestep, group] Julia matrices for clarity; serialised below
    # in wire order: all-lower rows, then all-median, then all-upper.
    lower_mat = Matrix{Float32}(undef, n_ts, N_GROUPS)
    median_mat = Matrix{Float32}(undef, n_ts, N_GROUPS)
    upper_mat = Matrix{Float32}(undef, n_ts, N_GROUPS)

    for g in 1:N_GROUPS
        for t in 1:n_ts
            vals = filter(!isnan, vec(results.group_cover[t, 1, g, :]))
            if isempty(vals)
                lower_mat[t, g] = NaN32
                median_mat[t, g] = NaN32
                upper_mat[t, g] = NaN32
            else
                q = quantile(vals, (0.025, 0.5, 0.975))
                lower_mat[t, g] = Float32(q[1])
                median_mat[t, g] = Float32(q[2])
                upper_mat[t, g] = Float32(q[3])
            end
        end
    end

    # ---- Serialise to wire layout ----
    buf = IOBuffer()
    buf.append = true  # ensure writes append

    # n_valid_runs (u32)
    write(buf, UInt32(n_valid))

    # covers[MAX_RUNS * N_TIMESTEPS] f32, run-major, zero-padded for unused runs
    covers_flat = zeros(Float32, MAX_RUNS * N_TIMESTEPS)
    for (col, r) in enumerate(valid_indices)
        base = (col - 1) * N_TIMESTEPS
        for t in 1:n_ts
            covers_flat[base + t] = Float32(results.cover[t, 1, r])
        end
    end
    write(buf, covers_flat)

    # summary: lower, median, upper — each [[f32; N_GROUPS]; N_TIMESTEPS] row-major
    # i.e. [t=1,g=1..5], [t=2,g=1..5], ..., [t=N_TIMESTEPS,g=1..5]
    for stat_mat in (lower_mat, median_mat, upper_mat)
        for t in 1:n_ts
            for g in 1:N_GROUPS
                write(buf, stat_mat[t, g])
            end
        end
    end

    # dhw[N_TIMESTEPS] f32 -- per-timestep DHW magnitude for the single simulated site
    for t in 1:n_ts
        write(buf, dhw_mat[t, 1])
    end

    result = take!(buf)
    length(result) == WORKER_RESULT_BYTES || error(
        "result size mismatch: wrote $(length(result)), expected $WORKER_RESULT_BYTES"
    )
    return result
end

# ---------------------------------------------------------------------------
# Zero-filled error result (n_valid_runs = 0) returned on simulation failure
# ---------------------------------------------------------------------------
function error_result()::Vector{UInt8}
    return zeros(UInt8, WORKER_RESULT_BYTES)
end

# ---------------------------------------------------------------------------
# Main worker loop
# ---------------------------------------------------------------------------
function run(args::Vector{String})::Cint
    if length(args) < 2
        println(
            Core.stderr,
            "[kora-worker] usage: kora-worker <growth_model_path> <survival_model_path>"
        )
        return Cint(1)
    end

    _ensure_blocking_stdin()

    # Load models once at startup
    try
        gm = Kora.load_models(args[1])::Kora.PolyGrowthModel{Float32}
        sm = Kora.load_models(args[2])::Kora.PolySurvivalModel{Float32}
        _growth_ref[] = gm
        _survival_ref[] = sm
        Kora._set_models!(gm, sm)
        println(Core.stderr, "[kora-worker] models loaded OK")
        flush(Core.stderr)
    catch
        println(Core.stderr, "[kora-worker] failed to load models")
        return Cint(1)
    end

    # Signal readiness to coordinator via stderr (stdout is binary-only)
    println(Core.stderr, "[kora-worker] READY")
    flush(Core.stderr)

    while true
        bytes = read_exact_stdin(WORKER_PARAMS_BYTES)
        bytes === nothing && break  # clean EOF — coordinator closed the pipe

        result = try
            p = parse_params(bytes)
            run_simulation(p)
        catch
            println(Core.stderr, "[kora-worker] simulation error")
            flush(Core.stderr)
            error_result()
        end

        write(Core.stdout, result)
        flush(Core.stdout)
    end

    println(Core.stderr, "[kora-worker] stdin closed, exiting")
    return Cint(0)
end

end  # module KoraWorker

# ---------------------------------------------------------------------------
# juliac --output-exe entry point (top-level, not inside a module)
# ---------------------------------------------------------------------------
function main(ARGS::Vector{String})::Cint
    return KoraWorker.run(ARGS)
end
Base.@main
