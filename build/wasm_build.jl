# build/wasm_build.jl -- Compile Kora.jl to WebAssembly via WasmTarget.jl.
#
# Usage (from Kora.jl root):
#   julia --project=. build/wasm_build.jl [output_dir]
#
# Prereqs:
#   - WasmTarget.jl added to the project:  julia --project=. -e 'import Pkg; Pkg.add(url="https://github.com/GroupTherapyOrg/WasmTarget.jl")'
#   - wasm-tools in PATH:                  cargo install wasm-tools
#
# Known blockers (see .claude/plans/wasm/):
#   - Distributions / LogNormal in ReefState.jl:496-501 not yet replaced
#   - Bootstrap, StatsBase, CSV/DataFrames, DiskArrays dep audit pending
#
# strict=false is intentional: emits traps for unsupported constructs instead of
# aborting, so the full blocker list surfaces in one pass rather than one at a time.

using WasmTarget
using Kora

# ---------------------------------------------------------------------------
# WASM bridge module
#
# Adapts the native C-pointer interface from bridge_aot.jl to types that
# WasmGC can represent:
#   - Ptr{UInt8} / unsafe_string  -> String
#   - Ptr{Float32} output buffers -> returned Vector{Float32} (flat layout)
#   - filesystem open()           -> embedded at compile time via __init__ using
#                                    load_models_from_string
# ---------------------------------------------------------------------------

module KoraWasm

using Kora
using Random

const _growth_ref = Ref{Union{Nothing,Kora.PolyGrowthModel{Float32}}}(nothing)
const _survival_ref = Ref{Union{Nothing,Kora.PolySurvivalModel{Float32}}}(nothing)
const _dhw_ref = Ref{Union{Nothing,Matrix{Float32}}}(nothing)
const _init_area_ref = Ref{Float32}(0.0f0)

const _deploy_vols_ref = Ref{NTuple{5,UInt32}}((
    UInt32(0), UInt32(0), UInt32(0), UInt32(0), UInt32(0)
))
const _deploy_start_ref = Ref{UInt32}(UInt32(1))
const _deploy_cadence_ref = Ref{UInt32}(UInt32(1))
const _depth_ref = Ref{Float32}(9.0f0)
const _deploy_dhw_tol_ref = Ref{Float32}(0.0f0)

# Stores the last kw_run_reef result so kw_result_get can serve it element-by-element.
# Avoids returning a WasmGC Vector ref to JS (which is opaque and unreadable directly).
const _last_result_ref = Ref{Union{Nothing,Vector{Float32}}}(nothing)

const _N_TIMESTEPS = Int32(75)
const _N_GROUPS = Int32(5)

# Stub: kw_load_models is kept as a stub (returns -2) so the caller can
# distinguish "not available in WASM build" from other errors.  The actual
# model data is embedded at compile time via __init__ below.
function kw_load_models(::String, ::String)::Int32
    return Int32(-2)
end

const _GROWTH_JSON = read(
    joinpath(@__DIR__, "../assets/models/offshore_north_growth_models.json"), String
)
const _SURVIVAL_JSON = read(
    joinpath(@__DIR__, "../assets/models/offshore_north_survival_models.json"), String
)

function __init__()
    gm = Kora.load_models_from_string(_GROWTH_JSON)::Kora.PolyGrowthModel{Float32}
    sm = Kora.load_models_from_string(_SURVIVAL_JSON)::Kora.PolySurvivalModel{Float32}
    _growth_ref[] = gm
    _survival_ref[] = sm
end

function kw_set_deployment(
    vol0::UInt32, vol1::UInt32, vol2::UInt32, vol3::UInt32, vol4::UInt32,
    start_year::UInt32,
    cadence_years::UInt32,
    depth_m::Float32,
    dhw_tol::Float32
)::Int32
    _deploy_vols_ref[] = (vol0, vol1, vol2, vol3, vol4)
    _deploy_start_ref[] = start_year
    _deploy_cadence_ref[] = cadence_years
    _depth_ref[] = depth_m
    _deploy_dhw_tol_ref[] = dhw_tol
    return Int32(0)
end

# kw_run_reef_d: deployment-aware variant that accepts deployment params directly.
# Replaces the kw_set_deployment + kw_run_reef two-call pattern in JS.
# kw_set_deployment conflicts with kw_run_reef during WasmTarget bundling (type-index
# renumbering bug), so we fold all deployment params into a single entry point instead.
function kw_run_reef_d(
    area_m2::Float32, init_cover_pct::Float32, n_runs::UInt32,
    vol0::UInt32, vol1::UInt32, vol2::UInt32, vol3::UInt32, vol4::UInt32,
    deploy_start::UInt32, deploy_cadence::UInt32,
    depth_m::Float32,
    dhw_tol::Float32
)::Int32
    _deploy_vols_ref[] = (vol0, vol1, vol2, vol3, vol4)
    _deploy_start_ref[] = deploy_start
    _deploy_cadence_ref[] = deploy_cadence
    # Inlined body of kw_run_reef — calling kw_run_reef directly triggers WasmTarget's
    # type-index renumbering bug when kw_run_reef is included as a private helper.
    n_ts = Int(_N_TIMESTEPS)
    n_groups = Int(_N_GROUPS)
    n_members = Int(n_runs)

    gm = _growth_ref[]
    gm === nothing && (_last_result_ref[]=nothing; return Int32(0))
    sm = _survival_ref[]
    sm === nothing && (_last_result_ref[]=nothing; return Int32(0))

    if _dhw_ref[] === nothing || _init_area_ref[] != area_m2
        _init_area_ref[] = area_m2
        _dhw_ref[] = Kora._wasm_generate_dhw(n_ts, 1)
    end
    dhw_mat = _dhw_ref[]::Matrix{Float32}

    reef = Kora.initialize_reef(n_ts, 1, Float64(area_m2), 10, Float64(depth_m), gm, sm)
    Kora._wasm_init_coral_pop!(reef)

    vols = _deploy_vols_ref[]
    start = Int(_deploy_start_ref[])
    cadence = Int(_deploy_cadence_ref[])
    if start >= 1 && cadence >= 1
        ts::Int = start
        while ts <= n_ts
            @inbounds for grp::Int in 1:5
                reef.deployment_times[ts, 1, grp] = Float32(vols[grp])
            end
            ts += cadence
        end
    end

    mean_cov = Float64(Kora.mean_colony_cover_m2())
    target_cover_m2 = (Float64(init_cover_pct) / 100.0) * Float64(area_m2)
    target_pop = max(5, ceil(Int64, target_cover_m2 / mean_cov))
    pop_density = Float64(target_pop) / Float64(area_m2)
    params = Matrix{Float64}(undef, 6, n_members)
    @inbounds for j::Int in 1:n_members
        params[1, j] = pop_density
    end
    @inbounds for i::Int in 2:6, j::Int in 1:n_members
        params[i, j] = 0.2
    end

    results = Kora._wasm_run_ensemble!(reef, dhw_mat, params, dhw_tol)

    n_valid::Int = 0
    @inbounds for r::Int in 1:n_members
        has_nan = false
        @inbounds for t2::Int in 1:n_ts
            if isnan(results.cover[t2, 1, r])
                has_nan = true
                break
            end
        end
        if !has_nan
            n_valid += 1
        end
    end

    buf = Vector{Float32}(undef, 2 + n_ts + n_ts * n_valid + n_ts * n_groups * 3)
    off = 1
    @inbounds buf[off] = Float32(n_ts);
    off += 1
    @inbounds buf[off] = Float32(n_valid);
    off += 1

    @inbounds for t::Int in 1:n_ts
        buf[off] = dhw_mat[t, 1];
        off += 1
    end

    @inbounds for r::Int in 1:n_members
        has_nan = false
        @inbounds for t2::Int in 1:n_ts
            if isnan(results.cover[t2, 1, r])
                has_nan = true
                break
            end
        end
        if !has_nan
            @inbounds for t::Int in 1:n_ts
                buf[off] = Float32(results.cover[t, 1, r]);
                off += 1
            end
        end
    end

    scratch = Vector{Float32}(undef, n_members)
    @inbounds for g::Int in 1:n_groups
        @inbounds for t::Int in 1:n_ts
            n_ok::Int = 0
            @inbounds for r::Int in 1:n_members
                v = results.group_cover[t, 1, g, r]
                if !isnan(v)
                    n_ok += 1
                    scratch[n_ok] = v
                end
            end
            lo, med, hi = if n_ok == 0
                NaN32, NaN32, NaN32
            else
                scratch_n = scratch[1:n_ok]
                sort!(scratch_n)
                _q(scratch_n, n_ok, 0.025f0),
                _q(scratch_n, n_ok, 0.5f0),
                _q(scratch_n, n_ok, 0.975f0)
            end
            @inbounds buf[off] = lo;
            off += 1
            @inbounds buf[off] = med;
            off += 1
            @inbounds buf[off] = hi;
            off += 1
        end
    end

    _last_result_ref[] = buf
    return Int32(length(buf))
end

# Linear-interpolation quantile on a pre-sorted slice scratch[1:n].
# Avoids Statistics.quantile, which pulls in iterator struct types that
# WasmTarget's compile_new cannot handle (BoundsError on 2-element field_types).
function _q(scratch::Vector{Float32}, n::Int, p::Float32)::Float32
    n == 1 && return scratch[1]
    h = p * Float32(n - 1)
    lo_i = floor(Int, h) + 1
    hi_i = min(lo_i + 1, n)
    frac = h - Float32(lo_i - 1)
    return scratch[lo_i] * (1.0f0 - frac) + scratch[hi_i] * frac
end

# Runs the reef simulation ensemble.
#
# Stores the result flat buffer in _last_result_ref[] and returns its length
# as Int32 (0 on error). JS reads elements via kw_result_get(i).
#
# Buffer layout (same as before, now accessed via kw_result_get):
#   [0]           n_ts    (Float32-encoded Int)
#   [1]           n_valid (Float32-encoded Int; valid ensemble members)
#   [2..2+n_ts)   DHW timeseries (single scenario column)
#   then          coral cover: n_valid * n_ts floats, run-major
#   then          stats: n_groups * n_ts * 3 floats [lower, median, upper] per group/timestep
function kw_run_reef(
    area_m2::Float32, init_cover_pct::Float32, n_runs::UInt32
)::Int32
    n_ts = Int(_N_TIMESTEPS)
    n_groups = Int(_N_GROUPS)
    n_members = Int(n_runs)

    gm = _growth_ref[]
    gm === nothing && (_last_result_ref[]=nothing; return Int32(0))
    sm = _survival_ref[]
    sm === nothing && (_last_result_ref[]=nothing; return Int32(0))

    if _dhw_ref[] === nothing || _init_area_ref[] != area_m2
        _init_area_ref[] = area_m2
        _dhw_ref[] = Kora._wasm_generate_dhw(n_ts, 1)
    end
    dhw_mat = _dhw_ref[]::Matrix{Float32}

    reef = Kora.initialize_reef(
        n_ts, 1, Float64(area_m2), 10, Float64(_depth_ref[]), gm, sm
    )
    Kora._wasm_init_coral_pop!(reef)

    vols = _deploy_vols_ref[]
    start = Int(_deploy_start_ref[])
    cadence = Int(_deploy_cadence_ref[])
    if start >= 1 && cadence >= 1
        ts::Int = start
        while ts <= n_ts
            @inbounds for grp::Int in 1:5
                reef.deployment_times[ts, 1, grp] = Float32(vols[grp])
            end
            ts += cadence
        end
    end

    dhw_tol = _deploy_dhw_tol_ref[]
    if dhw_tol > 0.0f0
        dhw_tol_sigma = 0.5f0
        start2 = Int(_deploy_start_ref[])
        cadence2 = Int(_deploy_cadence_ref[])
        for ts2 in start2:cadence2:n_ts, grp in 1:5
            if vols[grp] > 0
                reef.deployed_dhw_tolerances[ts2, 1, grp, 1] = dhw_tol
                reef.deployed_dhw_tolerances[ts2, 1, grp, 2] = dhw_tol_sigma
            end
        end
    end

    mean_cov = Float64(Kora.mean_colony_cover_m2())
    target_cover_m2 = (Float64(init_cover_pct) / 100.0) * Float64(area_m2)
    target_pop = max(5, ceil(Int64, target_cover_m2 / mean_cov))
    pop_density = Float64(target_pop) / Float64(area_m2)
    params = Matrix{Float64}(undef, 6, n_members)
    @inbounds for j::Int in 1:n_members
        params[1, j] = pop_density
    end
    @inbounds for i::Int in 2:6, j::Int in 1:n_members
        params[i, j] = 0.2
    end

    results = Kora._wasm_run_ensemble!(reef, dhw_mat, params, dhw_tol)

    # Pass 1: count valid ensemble members (no NaN in cover) without Vector{Int} —
    # WasmGC boxes Int64 as GC refs in arrays so storing i64 into Vector{Int} triggers
    # a browser validation error.
    n_valid::Int = 0
    @inbounds for r::Int in 1:n_members
        has_nan = false
        @inbounds for t2::Int in 1:n_ts
            if isnan(results.cover[t2, 1, r])
                has_nan = true
                break
            end
        end
        if !has_nan
            n_valid += 1
        end
    end

    buf = Vector{Float32}(undef, 2 + n_ts + n_ts * n_valid + n_ts * n_groups * 3)
    off = 1
    @inbounds buf[off] = Float32(n_ts);
    off += 1
    @inbounds buf[off] = Float32(n_valid);
    off += 1

    @inbounds for t::Int in 1:n_ts
        buf[off] = dhw_mat[t, 1];
        off += 1
    end

    # Pass 2: write cover for valid runs (re-check NaN rather than storing indices).
    @inbounds for r::Int in 1:n_members
        has_nan = false
        @inbounds for t2::Int in 1:n_ts
            if isnan(results.cover[t2, 1, r])
                has_nan = true
                break
            end
        end
        if !has_nan
            @inbounds for t::Int in 1:n_ts
                buf[off] = Float32(results.cover[t, 1, r]);
                off += 1
            end
        end
    end

    # Scratch buffer for per-timestep quantiles.
    # sort! on a copy (scratch[1:n_ok]) avoids SubArray — view() is WasmGC-incompatible.
    scratch = Vector{Float32}(undef, n_members)
    @inbounds for g::Int in 1:n_groups
        @inbounds for t::Int in 1:n_ts
            n_ok::Int = 0
            @inbounds for r::Int in 1:n_members
                v = results.group_cover[t, 1, g, r]
                if !isnan(v)
                    n_ok += 1
                    scratch[n_ok] = v
                end
            end
            lo, med, hi = if n_ok == 0
                NaN32, NaN32, NaN32
            else
                scratch_n = scratch[1:n_ok]
                sort!(scratch_n)
                _q(scratch_n, n_ok, 0.025f0),
                _q(scratch_n, n_ok, 0.5f0),
                _q(scratch_n, n_ok, 0.975f0)
            end
            @inbounds buf[off] = lo;
            off += 1
            @inbounds buf[off] = med;
            off += 1
            @inbounds buf[off] = hi;
            off += 1
        end
    end

    _last_result_ref[] = buf
    return Int32(length(buf))
end

# Returns the i-th float (0-indexed) from the last kw_run_reef result.
# Returns 0f0 if the index is out of range or no result is stored.
function kw_result_get(i::Int32)::Float32
    r = _last_result_ref[]
    r === nothing && return 0.0f0
    idx = Int(i) + 1
    (idx < 1 || idx > length(r)) && return 0.0f0
    return @inbounds r[idx]
end

# ---------------------------------------------------------------------------
# Sub-probes: isolate which Kora sub-call triggers WasmTarget BoundsError
#
# Each probe wraps one logical stage of kw_run_reef so the failing stage can
# be identified by running the probe list and observing which first fails.
# ---------------------------------------------------------------------------

# Stage 1: DHW generation only.
function probe_dhw(n_ts::Int32)::Int32
    dhw = Kora._wasm_generate_dhw(Int64(n_ts), Int64(1))
    return Int32(size(dhw, 1))
end

# Stage 0c: update_wild_sample! — stores a Vector{Float32} into wild_population.
function probe_update_wild(area::Float32)::Int32
    gm = _growth_ref[]
    gm === nothing && return Int32(-1)
    sm = _survival_ref[]
    sm === nothing && return Int32(-1)
    reef = Kora.initialize_reef(;
        n_timesteps=Int(_N_TIMESTEPS), n_locs=1, area=Float64(area), density=10,
        growth_models=gm, survival_models=sm
    )
    pop = Vector{Float32}([5.0f0, 10.0f0])
    Kora.update_wild_sample!(reef, 1, 1, 1, pop)
    return Int32(0)
end

# Stage 0d: DHW tolerance broadcast — 4D array slice assignment.
function probe_tolerances(area::Float32)::Int32
    gm = _growth_ref[]
    gm === nothing && return Int32(-1)
    sm = _survival_ref[]
    sm === nothing && return Int32(-1)
    reef = Kora.initialize_reef(;
        n_timesteps=Int(_N_TIMESTEPS), n_locs=1, area=Float64(area), density=10,
        growth_models=gm, survival_models=sm
    )
    # Explicit loops: view()/broadcast creates SubArray structs that WasmTarget's compile_new cannot handle
    for loc in 1:size(reef.wild_dhw_tolerances, 2)
        reef.wild_dhw_tolerances[1, loc, 1, 1] = 3.751612251f0
    end
    for ts in 1:size(reef.wild_dhw_tolerances, 1),
        loc in 1:size(reef.wild_dhw_tolerances, 2)

        reef.wild_dhw_tolerances[ts, loc, 1, 2] = 2.904433676f0
    end
    return Int32(0)
end

# Stage 0a: bin_edges() only — tests size_classes.jl with no RNG or reef state.
function probe_bin_edges()::Int32
    edges = Kora.bin_edges()
    return Int32(size(edges, 1))
end

# Stage 0b: _sample_lognormal_bounded only — tests the randn-based sampling loop.
function probe_sample(n::Int32)::Int32
    rng = Random.default_rng()
    out = Kora._sample_lognormal_bounded(2.238f0, 0.749f0, 0.0, 160.0, Int(n), rng)
    return Int32(length(out))
end

# Stage 2a: initialize_reef only (no coral population).
function probe_init_reef(area::Float32)::Int32
    gm = _growth_ref[]
    gm === nothing && return Int32(-1)
    sm = _survival_ref[]
    sm === nothing && return Int32(-1)
    _ = Kora.initialize_reef(;
        n_timesteps=Int(_N_TIMESTEPS), n_locs=1, area=Float64(area), density=10,
        growth_models=gm, survival_models=sm
    )
    return Int32(0)
end

# Stage 2b: initialize_reef + initialize_coral_population! (no ensemble run).
function probe_init(area::Float32)::Int32
    gm = _growth_ref[]
    gm === nothing && return Int32(-1)
    sm = _survival_ref[]
    sm === nothing && return Int32(-1)
    reef = Kora.initialize_reef(;
        n_timesteps=Int(_N_TIMESTEPS), n_locs=1, area=Float64(area), density=10,
        growth_models=gm, survival_models=sm
    )
    Kora.initialize_coral_population!(reef)
    return Int32(0)
end

# Stage 3: full ensemble run, no post-processing (no quantile / buffer writes).
function probe_run(area::Float32, n::UInt32)::Int32
    gm = _growth_ref[]
    gm === nothing && return Int32(-1)
    sm = _survival_ref[]
    sm === nothing && return Int32(-1)
    dhw = Kora._wasm_generate_dhw(Int64(_N_TIMESTEPS), Int64(1))
    reef = Kora.initialize_reef(;
        n_timesteps=Int(_N_TIMESTEPS), n_locs=1, area=Float64(area), density=10,
        growth_models=gm, survival_models=sm
    )
    Kora.initialize_coral_population!(reef)
    n_members = Int(n)
    params = Matrix{Float64}(undef, 6, n_members)
    for j::Int in 1:n_members
        params[1, j] = 1.0
    end
    for i::Int in 2:6, j::Int in 1:n_members
        params[i, j] = 0.2
    end
    _ = Kora._wasm_run_ensemble!(reef, dhw, params, 0.0f0)
    return Int32(0)
end

end  # module KoraWasm

# ---------------------------------------------------------------------------
# Compile
#
# Each entry point is probed individually first so failures are attributed to
# a specific function rather than silently aborting the whole build.
# Successfully compiled entries are then bundled into a single .wasm module.
# ---------------------------------------------------------------------------

out_dir = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "dist")
out_path = joinpath(out_dir, "kora.wasm")
mkpath(out_dir)

all_entries = [
    ("kw_load_models", KoraWasm.kw_load_models, (String, String)),
    (
        "kw_set_deployment",
        KoraWasm.kw_set_deployment,
        (UInt32, UInt32, UInt32, UInt32, UInt32, UInt32, UInt32)
    ),
    ("kw_run_reef", KoraWasm.kw_run_reef, (Float32, Float32, UInt32)),
    (
        "kw_run_reef_d",
        KoraWasm.kw_run_reef_d,
        (Float32, Float32, UInt32, UInt32, UInt32, UInt32, UInt32, UInt32, UInt32, UInt32)
    ),
    ("kw_result_get", KoraWasm.kw_result_get, (Int32,)),
    # Sub-probes: compiled but excluded from the final bundle (see below).
    # Run these to isolate which Kora stage triggers the WasmTarget BoundsError.
    ("probe_bin_edges", KoraWasm.probe_bin_edges, ()),
    ("probe_sample", KoraWasm.probe_sample, (Int32,)),
    ("probe_update_wild", KoraWasm.probe_update_wild, (Float32,)),
    ("probe_tolerances", KoraWasm.probe_tolerances, (Float32,)),
    ("probe_dhw", KoraWasm.probe_dhw, (Int32,)),
    ("probe_init_reef", KoraWasm.probe_init_reef, (Float32,)),
    ("probe_init", KoraWasm.probe_init, (Float32,)),
    ("probe_run", KoraWasm.probe_run, (Float32, UInt32))
]

# Sub-probes are diagnostic only — never included in the shipped .wasm bundle.
# kw_load_models is also excluded: it is a permanent stub (returns -2) and its
# String arguments cause a WasmTarget type-registry index mismatch when bundled
# with the numeric entry points. It is provided by the JS host instead.
const _PROBE_NAMES = Set(["probe_bin_edges", "probe_sample", "probe_update_wild",
    "probe_tolerances",
    "probe_dhw", "probe_init_reef", "probe_init", "probe_run",
    "kw_load_models",
    "kw_set_deployment",  # conflicts with kw_run_reef during bundling (type-index renumbering)
    "kw_run_reef"])        # superseded by kw_run_reef_d which folds deployment params in

good_entries = Tuple{Any,Any}[]
passed_probes = Set{String}()  # probes that compiled — skip these in diagnostics

for (name, fn, types) in all_entries
    @info "Probing $name ..."
    try
        WasmTarget.compile_multi(
            [(fn, types)]; strict=false, validate=false, discovery=:trim
        )
        if name in _PROBE_NAMES
            push!(passed_probes, name)
        else
            push!(good_entries, (fn, types))
        end
        @info "  $name: OK"
    catch e
        bt = catch_backtrace()
        @warn "  $name: FAILED\n$(sprint(showerror, e, bt))"
    end
end

if isempty(good_entries)
    error("All entry points failed to compile — nothing to emit.")
end

# Diagnostic: walk the full typed call graph of failing entry points and report
# any :new expression whose IR argument count != Julia fieldcount.  That mismatch
# is the WasmTarget BoundsError trigger (statements.jl:2004: field_types[i] OOB).
#
# Uses Base.code_typed_by_type(mi.specTypes) to follow :invoke edges via the
# MethodInstance directly — avoids the getfield(module, name) reconstruction that
# fails for closures and keyword-argument wrapper methods.
#
# Output format:
#   [depth] :new(OffendingType) — IR n_ir args vs Julia n_jl fields

const _MAX_DIAG_DEPTH = 12

function _scan_new_mismatches(fn, types; max_depth=_MAX_DIAG_DEPTH)
    visited = Set{Any}()   # MethodInstance or :toplevel sentinel
    found = Ref(false)

    function check_ci(ci, depth)
        depth > max_depth && return nothing
        for (i, stmt) in enumerate(ci.code)
            if stmt isa Expr && stmt.head === :new
                T = stmt.args[1]
                if T isa DataType && !isabstracttype(T) && !(T <: Array)
                    n_ir = length(stmt.args) - 1
                    n_jl = try
                        fieldcount(T)
                    catch
                        ;
                        -1
                    end
                    if n_ir != n_jl
                        indent = "  " ^ depth
                        @warn "$(indent)[depth=$depth] :new($T) — IR $n_ir args, Julia $n_jl fields  (stmt $i)"
                        found[] = true
                    end
                end
            end

            # Follow :invoke edges via the MethodInstance directly.
            if stmt isa Expr && stmt.head === :invoke
                mi = stmt.args[1]
                if mi isa Core.MethodInstance && !(mi in visited)
                    push!(visited, mi)
                    try
                        for (sub_ci, _) in
                            Base.code_typed_by_type(mi.specTypes; optimize=true)
                            check_ci(sub_ci, depth + 1)
                        end
                    catch
                    end
                end
            end
        end
    end

    push!(visited, :toplevel)
    try
        for (ci, _) in code_typed(fn, types; optimize=true)
            check_ci(ci, 1)
        end
    catch e
        @warn "  code_typed failed: $(sprint(showerror, e))"
    end
    found[] ||
        @info "  no :new mismatches found in call graph (visited $(length(visited)) methods)"
end

for (name, fn, types) in all_entries
    # Skip entries that compiled successfully (probes or real entries).
    (name in passed_probes || any(e -> e === (fn, types), good_entries)) && continue
    @info "Diagnosing :new call-graph for $name ..."
    try
        _scan_new_mismatches(fn, types)
    catch e2
        @warn "  diagnostic failed: $(sprint(showerror, e2))"
    end
end

@info "Bundling $(length(good_entries))/$(length(all_entries)) entry points ..."
t = @elapsed begin
    bytes = WasmTarget.compile_multi(
        good_entries; strict=false, validate=false, discovery=:trim)
end
@info "Done in $(round(t; digits=1))s -- $(length(bytes) >>> 10) KB"
write(out_path, bytes)
@info "Written -> $out_path"
