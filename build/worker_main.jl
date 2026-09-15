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
# WorkerSimParams wire layout (little-endian, 1676 bytes, no padding):
#   reef_area_m2:            f32  offset  0
#   init_cover_pct:          f32  offset  4
#   deploy_volumes[5]:       u32  offset  8  (20 bytes)
#   deploy_start_year:       u32  offset 28
#   deploy_cadence_years:    u32  offset 32
#   depth_m:                 u32  offset 36
#   deploy_dhw_tolerance:    f32  offset 40
#   dhw_seed:                u32  offset 44
#   init_group_fraction[5]:  f32  offset 48  (20 bytes)
#   dhw_override[MAX_TIMESTEPS]: f32 offset 68 (1200 bytes)   # Part 3, §6.2 option (b)
#   dhw_override_active:     u32  offset 1268 (4 bytes)
#   init_size_class_fraction[35]: f32 offset 1272 (140 bytes)   # Part 5 v2
#   init_size_class_active:  u32  offset 1412 (4 bytes)
#   n_timesteps:             u32  offset 1416 (4 bytes)   # Part 10 -- requested run length
#   Total: 1420 bytes
#
# NOTE: WorkerSimParams is a superset of sim-types/src/wire.rs WireSimParams —
# it adds depth_m and deploy_dhw_tolerance. wire.rs's WireSimParams matches
# this layout as of kora-app's feat/web-backend branch (Part 10: both now
# size their time-indexed fields to MAX_TIMESTEPS and carry an explicit
# n_timesteps -- see kora-app's `.claude/plans/revise-ui.md` §10, OQ-7 option (a)).
#
# WorkerEnsembleResult wire layout (little-endian, no padding):
#   n_valid_runs:                        u32  offset       0  (4 bytes)
#   covers[MAX_RUNS * MAX_TIMESTEPS]:    f32  offset       4  (120000 bytes), run-major
#   summary.lower[MAX_TIMESTEPS][N_GROUPS]:  f32  offset  120004  (6000 bytes)
#   summary.median[MAX_TIMESTEPS][N_GROUPS]: f32  offset  126004  (6000 bytes)
#   summary.upper[MAX_TIMESTEPS][N_GROUPS]:  f32  offset  132004  (6000 bytes)
#   dhw[MAX_TIMESTEPS]:                   f32  offset  138004  (1200 bytes)
#   tolerance.lower[MAX_TIMESTEPS][N_GROUPS]:  f32  offset  139204  (6000 bytes)
#   tolerance.median[MAX_TIMESTEPS][N_GROUPS]: f32  offset  145204  (6000 bytes)
#   tolerance.upper[MAX_TIMESTEPS][N_GROUPS]:  f32  offset  151204  (6000 bytes)
#   n_timesteps:                          u32  offset  157204  (4 bytes)   # Part 10 -- actual run length
#   Total: 157208 bytes
#
# summary / tolerance layout mirrors WireGroupSummary in wire.rs:
#   [[f32; N_GROUPS]; MAX_TIMESTEPS] = row-major with timestep as outer index.
#   Only the first n_timesteps rows/entries of every time-indexed field above
#   are meaningful -- the rest is zero padding out to MAX_TIMESTEPS capacity.
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
const N_SIZES = 7  # Part 5 v2 -- must stay in sync with sim-types/src/results.rs N_SIZES
# Default run length, used only as a defensive fallback when a request's
# n_timesteps field is missing/zero (shouldn't happen from a current client).
const N_TIMESTEPS = 75
# Capacity every time-indexed wire field is sized to (Part 10). Must match
# sim-types/src/results.rs::MAX_TIMESTEPS / wire.rs::WIRE_MAX_TIMESTEPS.
const MAX_TIMESTEPS = 300
const MAX_RUNS = 100

# 4 scalar fields (2x f32 + 5x u32 deploy_volumes + u32 start + u32 cadence) plus
# depth_m (u32), dhw_tol (f32), dhw_seed (u32), init_group_fraction (5x f32),
# dhw_override (MAX_TIMESTEPS x f32), dhw_override_active (u32),
# init_size_class_fraction (N_GROUPS*N_SIZES x f32), init_size_class_active (u32),
# n_timesteps (u32).
const WORKER_PARAMS_BYTES = 4 + 4 + N_GROUPS * 4 + 4 + 4 + 4 + 4 + 4 + N_GROUPS * 4 + MAX_TIMESTEPS * 4 + 4 +
                             N_GROUPS * N_SIZES * 4 + 4 + 4   # = 1420

# u32 n_valid + [MAX_RUNS * MAX_TIMESTEPS] f32 covers + [MAX_TIMESTEPS * N_GROUPS * 3] f32 summary
# + [MAX_TIMESTEPS] f32 dhw + [MAX_TIMESTEPS * N_GROUPS * 3] f32 tolerance + u32 n_timesteps
const WORKER_RESULT_BYTES = 4 + MAX_RUNS * MAX_TIMESTEPS * 4 + N_GROUPS * MAX_TIMESTEPS * 3 * 4 + MAX_TIMESTEPS * 4 + N_GROUPS * MAX_TIMESTEPS * 3 * 4 + 4  # = 157208

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
# from kora-server) sets pipe stdio to O_NONBLOCK on BOTH ends -- stdin (fd
# 0) and stdout (fd 1). A blocking-style raw ccall(:read)/ccall(:write) would
# then intermittently see EAGAIN (errno -11): on fd 0, indistinguishable from
# EOF (ret <= 0) unless we either check errno or just clear O_NONBLOCK once
# up front; on fd 1, a write() larger than the pipe's buffer capacity (65536
# bytes on Linux -- WORKER_RESULT_BYTES is 157208) does a partial write and
# then EAGAINs on the retry for the remainder instead of blocking until the
# reader drains it, which write_exact_stdout's retry loop can't tell apart
# from a real fatal write error. Clearing O_NONBLOCK once up front keeps both
# read_exact_stdin's "ret <= 0 means EOF" and write_exact_stdout's "ret <= 0
# means fatal" logic correct for both invocation styles (`< file` / `> file`
# redirection, where the fd was already blocking, and a live coordinator
# pipe).
function _ensure_blocking_fd(fd::Cint)::Nothing
    F_GETFL = Cint(3)
    F_SETFL = Cint(4)
    O_NONBLOCK = Cint(0o4000)
    flags = ccall(:fcntl, Cint, (Cint, Cint), fd, F_GETFL)
    flags >= 0 || return nothing
    if (flags & O_NONBLOCK) != 0
        ccall(:fcntl, Cint, (Cint, Cint, Cint), fd, F_SETFL, flags & ~O_NONBLOCK)
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

# Mirrors read_exact_stdin: a raw write(2) to a pipe is only guaranteed to
# accept up to the pipe's buffer capacity (65536 bytes on Linux) in one call
# -- for a payload larger than that (WORKER_RESULT_BYTES is 157208), the
# kernel does a partial write and the caller must retry with the remainder.
# `write(Core.stdout, result)` doesn't loop for this: Core.stdout resolves to
# a bare-bones IO under --trim=safe (same reason read_exact_stdin can't use
# ordinary `read(stdin, ...)`), and a single un-retried write silently
# truncates output at the pipe capacity boundary. The caller then loops back
# to read_exact_stdin waiting for the next request, so the truncation isn't
# even visible here -- it surfaces only as the client hanging forever short
# of WORKER_RESULT_BYTES.
function write_exact_stdout(bytes::Vector{UInt8})::Nothing
    n = length(bytes)
    total = 0
    while total < n
        ret = GC.@preserve bytes ccall(
            :write, Cssize_t,
            (Cint, Ptr{UInt8}, Csize_t),
            Cint(1), pointer(bytes, total + 1), Csize_t(n - total)
        )
        ret <= Cssize_t(0) && error("write to stdout failed")
        total += Int(ret)
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Parse WorkerSimParams from 516 raw bytes (little-endian field order)
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
    init_group_fraction = ntuple(_ -> read(io, Float32), N_GROUPS)
    # `Val(...)` (not a plain Integer) is required above 10 elements -- Base's
    # `ntuple(f, n::Integer)` only stays inferrable as a concrete `NTuple` for
    # n<=10 (an unrolled fast path); past that it returns the abstract
    # `Tuple{Vararg{T}}`, which juliac's --trim=safe verifier can't resolve a
    # downstream call against (caught by the AOT build, not by `julia --check-bounds`
    # or the test suite -- neither exercises trim verification).
    dhw_override = ntuple(_ -> read(io, Float32), Val(MAX_TIMESTEPS))
    dhw_override_active = read(io, UInt32)
    init_size_class_fraction = ntuple(_ -> read(io, Float32), Val(N_GROUPS * N_SIZES))
    init_size_class_active = read(io, UInt32)
    n_timesteps = read(io, UInt32)
    return (;
        reef_area_m2,
        init_cover_pct,
        deploy_volumes,
        deploy_start_year,
        deploy_cadence_years,
        depth_m,
        deploy_dhw_tolerance,
        dhw_seed,
        init_group_fraction,
        dhw_override,
        dhw_override_active,
        init_size_class_fraction,
        init_size_class_active,
        n_timesteps
    )
end

# ---------------------------------------------------------------------------
# Simulation helpers (adapted from bridge_aot.jl)
# ---------------------------------------------------------------------------
function _build_ensemble_params(
    area_m2::Float32, init_cover_pct::Float32, n_members::Int,
    group_fraction::NTuple{5,Float32},
    size_class_fraction::NTuple{35,Float32}, size_class_active::Bool
)::Matrix{Float64}
    mean_cov = Float64(Kora.mean_colony_cover_m2())
    target_cover_m2 = (Float64(init_cover_pct) / 100.0) * Float64(area_m2)
    target_pop = max(5, ceil(Int64, target_cover_m2 / mean_cov))
    pop_density = Float64(target_pop) / Float64(area_m2)

    n_rows = size_class_active ? 6 + N_GROUPS * N_SIZES : 6
    params = Matrix{Float64}(undef, n_rows, n_members)
    params[1, :] .= pop_density
    fr = collect(Float64.(group_fraction))
    s = sum(fr)
    fr = (s > 0 && isfinite(s)) ? fr ./ s : fill(0.2, 5)
    for g in 1:5
        params[1 + g, :] .= fr[g]
    end
    if size_class_active
        sc = collect(Float64.(size_class_fraction))
        for i in 1:(N_GROUPS * N_SIZES)
            params[6 + i, :] .= sc[i]
        end
    end
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
    # Part 10: run length now comes from the request, not a fixed constant.
    # A zero/missing n_timesteps (shouldn't happen from a current client)
    # falls back to the historical default rather than erroring.
    n_ts = p.n_timesteps == 0 ? N_TIMESTEPS : clamp(Int(p.n_timesteps), 1, MAX_TIMESTEPS)

    gm = _growth_ref[]
    sm = _survival_ref[]
    (gm === nothing || sm === nothing) && error("models not loaded")

    # Reuse cached DHW unless n_ts or the requested seed changed (same policy
    # as bridge_aot.jl, plus seed-awareness -- kora-server has no separate
    # regenerate-DHW endpoint, so a client asking for a different seed on an
    # otherwise ordinary /api/run_reef call is how "New DHW trajectory" is
    # expressed here).
    # Part 3, §6.2 option (b): a custom DHW trajectory rides inline in every
    # request. When active, use the first n_ts entries verbatim and skip the
    # seed cache entirely.
    if p.dhw_override_active != 0
        # Runtime-range slicing a fixed-size tuple (`p.dhw_override[1:n_ts]`)
        # can't stay concretely typed -- n_ts isn't known until runtime, so
        # juliac's --trim=safe verifier rejects it. Element-at-a-time Int
        # indexing into the tuple stays concrete regardless of n_ts.
        dhw_vec = Vector{Float32}(undef, n_ts)
        for t in 1:n_ts
            dhw_vec[t] = p.dhw_override[t]
        end
        dhw_mat = reshape(dhw_vec, n_ts, 1)
    else
        if _dhw_ref[] === nothing || _init_n_ts_ref[] != n_ts || _dhw_seed_ref[] != p.dhw_seed
            _init_n_ts_ref[] = n_ts
            _dhw_seed_ref[] = p.dhw_seed
            _dhw_ref[] = Kora.generate_example_dhw(n_ts, 1; rng=Xoshiro(Int(p.dhw_seed)))
        end
        dhw_mat = _dhw_ref[]::Matrix{Float32}
    end

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
    ensemble_params = _build_ensemble_params(
        p.reef_area_m2, p.init_cover_pct, n_members, p.init_group_fraction,
        p.init_size_class_fraction, p.init_size_class_active != 0
    )
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

    # Wild-population mean DHW tolerance per group per timestep, same percentile
    # pass over ensemble members. results.wild_dhw_tolerances is
    # (n_ts, n_locs, n_groups, 2, n_members) with dim 4 = [mean, std]; take mean.
    tol_lower_mat = Matrix{Float32}(undef, n_ts, N_GROUPS)
    tol_median_mat = Matrix{Float32}(undef, n_ts, N_GROUPS)
    tol_upper_mat = Matrix{Float32}(undef, n_ts, N_GROUPS)

    for g in 1:N_GROUPS
        for t in 1:n_ts
            vals = filter(!isnan, vec(results.wild_dhw_tolerances[t, 1, g, 1, :]))
            if isempty(vals)
                tol_lower_mat[t, g] = NaN32
                tol_median_mat[t, g] = NaN32
                tol_upper_mat[t, g] = NaN32
            else
                q = quantile(vals, (0.025, 0.5, 0.975))
                tol_lower_mat[t, g] = Float32(q[1])
                tol_median_mat[t, g] = Float32(q[2])
                tol_upper_mat[t, g] = Float32(q[3])
            end
        end
    end

    # ---- Serialise to wire layout ----
    buf = IOBuffer()
    buf.append = true  # ensure writes append

    # n_valid_runs (u32)
    write(buf, UInt32(n_valid))

    # covers[MAX_RUNS * MAX_TIMESTEPS] f32, run-major, zero-padded for unused
    # runs AND for timesteps beyond n_ts (Part 10 capacity buffer).
    covers_flat = zeros(Float32, MAX_RUNS * MAX_TIMESTEPS)
    for (col, r) in enumerate(valid_indices)
        base = (col - 1) * MAX_TIMESTEPS
        for t in 1:n_ts
            covers_flat[base + t] = Float32(results.cover[t, 1, r])
        end
    end
    write(buf, covers_flat)

    # summary: lower, median, upper — each [[f32; N_GROUPS]; MAX_TIMESTEPS]
    # row-major, i.e. [t=1,g=1..5], [t=2,g=1..5], ..., [t=MAX_TIMESTEPS,g=1..5].
    # Rows beyond n_ts are zero padding out to capacity.
    for stat_mat in (lower_mat, median_mat, upper_mat)
        for t in 1:MAX_TIMESTEPS
            for g in 1:N_GROUPS
                write(buf, t <= n_ts ? stat_mat[t, g] : 0f0)
            end
        end
    end

    # dhw[MAX_TIMESTEPS] f32 -- per-timestep DHW magnitude for the single
    # simulated site, zero-padded beyond n_ts.
    for t in 1:MAX_TIMESTEPS
        write(buf, t <= n_ts ? dhw_mat[t, 1] : 0f0)
    end

    # tolerance: lower, median, upper — same [[f32; N_GROUPS]; MAX_TIMESTEPS]
    # row-major layout as summary, appended before n_timesteps so existing
    # offsets don't move.
    for stat_mat in (tol_lower_mat, tol_median_mat, tol_upper_mat)
        for t in 1:MAX_TIMESTEPS
            for g in 1:N_GROUPS
                write(buf, t <= n_ts ? stat_mat[t, g] : 0f0)
            end
        end
    end

    # n_timesteps (u32) -- actual run length, appended last so existing
    # offsets stay stable (Part 10).
    write(buf, UInt32(n_ts))

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

    _ensure_blocking_fd(Cint(0))
    _ensure_blocking_fd(Cint(1))

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

        write_exact_stdout(result)
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
