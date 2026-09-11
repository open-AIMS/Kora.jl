# bridge_aot.jl — Julia AOT bridge for Kora reef simulation
#
# C API signatures (for Rust FFI declarations):
#   int32_t kf_load_models(const uint8_t* growth_path, const uint8_t* surv_path);
#   int32_t kf_set_deployment(uint32_t vol0, uint32_t vol1, uint32_t vol2,
#                             uint32_t vol3, uint32_t vol4,
#                             uint32_t start_year, uint32_t cadence_years,
#                             float depth_m, float deploy_dhw_tolerance);
#   int32_t kf_new_dhw_trajectory();  // invalidate cached DHW; next kf_run_reef
#                                     // call regenerates a fresh climate sequence
#   int32_t kf_set_initial_cover(const float* group_fraction, int32_t n /* must be 5 */);
#   int32_t kf_set_dhw_trajectory(const float* values, int32_t n /* n<=0 clears */);
#   int32_t kf_run_reef(float area_m2, float init_cover_pct, uint32_t n_runs,
#                       uint32_t dhw_seed,
#                       float* dhw_out, int32_t dhw_cap,
#                       float* covers_out, int32_t covers_cap,
#                       float* lower_out, float* median_out, float* upper_out,
#                       float* tol_lower_out, float* tol_median_out, float* tol_upper_out,
#                       int32_t stats_cap,   // shared by the group-cover and tolerance out buffers
#                       int64_t* n_ts_out, int64_t* n_valid_out);

module KoraBridge

using Kora
using Random: Xoshiro
using Statistics: quantile

# Raw fd-2 write, bypassing Julia's IO system entirely -- needed because this
# runs from inside a trimmed/@ccallable exception handler where the ordinary
# runtime (stderr, println) may not be safe to call into.
#
# The underlying libc symbol differs by platform: Windows' CRT (io.h) exposes
# it as `_write(int fd, const void *buffer, unsigned int count)`; POSIX libc
# (glibc/musl) exposes the unprefixed `write(int fd, const void *buf, size_t
# count)` -- `_write` does not exist there and a ccall to it aborts the
# process ("could not load symbol"), taking down the *original* exception
# report with it. `@static if` resolves at parse/lowering time on whichever
# host is doing the AOT compile (build.sh vs build.ps1), leaving a single
# literal-symbol `ccall` in the compiled output -- `ccall`'s function-spec
# argument must be a literal, so this can't be a runtime-selected variable.
@static if Sys.iswindows()
    @inline function _raw_write(fd::Cint, buf, n::Integer)::Cint
        ccall(:_write, Cint, (Cint, Ptr{UInt8}, Cuint), fd, buf, Cuint(n))
    end
else
    @inline function _raw_write(fd::Cint, buf, n::Integer)::Cint
        ccall(:write, Cint, (Cint, Ptr{UInt8}, Csize_t), fd, buf, Csize_t(n))
    end
end

macro _write_stderr(msg)
    n = ncodeunits(msg)
    :(_raw_write(Int32(2), $msg, $n))
end

# Staging Ref: the exception is stored here before calling _kf_write_inner_exc,
# so the @ccallable helper takes no Any argument (which @ccallable forbids).
const _exc_stage = Ref{Any}(nothing)

# Store the exception in _exc_stage, print its type name, then call the
# @ccallable helper (no Any args → no verifier error).
# setindex! on Ref{Any} dispatches on the concrete first-arg type — fine.
macro _write_exception(e)
    quote
        let _e = $(esc(e))
            _exc_stage[] = _e
            _tcstr = ccall(:jl_typeof_str, Ptr{UInt8}, (Any,), _e)
            _tlen = ccall(:strlen, Csize_t, (Ptr{UInt8},), _tcstr)
            _raw_write(Int32(2), _tcstr, _tlen)
            # Direct Julia call (no Any args) — no verifier error, no Windows symbol
            # lookup issue.  The exception was already stored in _exc_stage above.
            _kf_write_inner_exc()
        end
    end
end

# Plain Julia function (not @ccallable — no need to export it to C).
# Zero-argument call from the macro means dispatch is unambiguous; verifier accepts it.
# Reads from _exc_stage and drills into CompositeException without jl_arrayref
# (not exported on Windows) by using Julia field + array indexing on narrowed types.
function _kf_write_inner_exc()::Nothing
    exc = _exc_stage[]
    if exc isa MethodError
        fn = (exc::MethodError).f
        fstr = ccall(:jl_typeof_str, Ptr{UInt8}, (Any,), fn)
        flen = ccall(:strlen, Csize_t, (Ptr{UInt8},), fstr)
        _raw_write(Int32(2), " on ", 4)
        _raw_write(Int32(2), fstr, flen)
    elseif exc isa CompositeException
        excs = (exc::CompositeException).exceptions
        if !isempty(excs)
            tfe = excs[1]
            if tfe isa TaskFailedException
                inner = (tfe::TaskFailedException).task.result
                icstr = ccall(:jl_typeof_str, Ptr{UInt8}, (Any,), inner)
                ilen = ccall(:strlen, Csize_t, (Ptr{UInt8},), icstr)
                _raw_write(Int32(2), " [inner: ", 9)
                _raw_write(Int32(2), icstr, ilen)
                if inner isa MethodError
                    fn = (inner::MethodError).f
                    fstr = ccall(:jl_typeof_str, Ptr{UInt8}, (Any,), fn)
                    flen = ccall(:strlen, Csize_t, (Ptr{UInt8},), fstr)
                    _raw_write(Int32(2), " on ", 4)
                    _raw_write(Int32(2), fstr, flen)
                end
                _raw_write(Int32(2), "]", 1)
            end
        end
    end
    return nothing
end

const _N_GROUPS = Int32(5)

const _growth_ref = Ref{Union{Nothing,Kora.PolyGrowthModel{Float32}}}(nothing)
const _survival_ref = Ref{Union{Nothing,Kora.PolySurvivalModel{Float32}}}(nothing)

# Cache DHW so every run batch uses the same climate forcing. Area does NOT
# invalidate this cache — DHW generation doesn't depend on area, and reef
# area changes must not silently change the climate trajectory. Regenerated
# only when timestep count or seed changes, models are (re)loaded, or a new
# trajectory is explicitly requested via kf_new_dhw_trajectory.
const _dhw_ref = Ref{Union{Nothing,Matrix{Float32}}}(nothing)
const _init_n_ts_ref = Ref{Int}(0)
const _dhw_seed_ref = Ref{UInt32}(0)

# Custom DHW/climate trajectory override (Part 3). When non-nothing and its
# length matches n_ts, kf_run_reef uses it verbatim and skips the seed cache /
# generate_example_dhw entirely. Set via kf_set_dhw_trajectory; cleared by
# kf_set_dhw_trajectory(n<=0), kf_new_dhw_trajectory, and kf_load_models.
const _dhw_override_ref = Ref{Union{Nothing,Vector{Float32}}}(nothing)

# Deployment schedule — set via kf_set_deployment before kf_run_reef.
# NTuple{5,UInt32}: corals/year for each of the 5 functional groups.
const _deploy_vols_ref = Ref{NTuple{5,UInt32}}((
    UInt32(0), UInt32(0), UInt32(0), UInt32(0), UInt32(0)
))
const _deploy_start_ref = Ref{UInt32}(UInt32(1))
const _deploy_cadence_ref = Ref{UInt32}(UInt32(1))
const _depth_ref = Ref{Float32}(9.0f0)
const _deploy_dhw_tol_ref = Ref{Float32}(0.0f0)

# Per-group initial-cover composition — set via kf_set_initial_cover before
# kf_run_reef. NTuple{5,Float32}: composition fraction for each of the 5
# functional groups (need not sum to 1; normalized in _build_ensemble_params).
const _group_fraction_ref = Ref{NTuple{5,Float32}}((0.2f0, 0.2f0, 0.2f0, 0.2f0, 0.2f0))

# Build ensemble params where all members share the same initial conditions
# (equal group proportions, cover-derived density) so CI-band spread at t=0
# reflects stochastic dynamics only, not variation in initial setup.
function _build_ensemble_params(
    area_m2::Float32, init_cover_pct::Float32, n_members::Int
)::Matrix{Float64}
    mean_cov = Float64(Kora.mean_colony_cover_m2())
    target_cover_m2 = (Float64(init_cover_pct) / 100.0) * Float64(area_m2)
    target_pop = max(5, ceil(Int64, target_cover_m2 / mean_cov))
    pop_density = Float64(target_pop) / Float64(area_m2)
    params = Matrix{Float64}(undef, 6, n_members)
    params[1, :] .= pop_density
    fr = collect(Float64.(_group_fraction_ref[]))
    s = sum(fr)
    fr = (s > 0 && isfinite(s)) ? fr ./ s : fill(0.2, 5)
    for g in 1:5
        params[1 + g, :] .= fr[g]
    end
    return params
end

# Populate reef_state.deployment_times from the global deployment schedule.
# All args are concrete scalars — no pointer reads, no dispatch issues.
function _apply_deployment!(reef::Kora.ReefState, n_ts::Int)::Nothing
    vols = _deploy_vols_ref[]
    start = Int(_deploy_start_ref[])
    cadence = Int(_deploy_cadence_ref[])
    if start >= 1 && cadence >= 1
        for ts::Int in start:cadence:n_ts
            for grp::Int in 1:5
                reef.deployment_times[ts, 1, grp] = Float32(vols[grp])
            end
        end
    end
    return nothing
end

Base.@ccallable function kf_load_models(
    growth_path::Ptr{UInt8},
    surv_path::Ptr{UInt8}
)::Int32
    try
        gp = unsafe_string(growth_path)
        sp = unsafe_string(surv_path)
        gm = Kora.load_models(gp)::Kora.PolyGrowthModel{Float32}
        sm = Kora.load_models(sp)::Kora.PolySurvivalModel{Float32}
        _growth_ref[] = gm
        _survival_ref[] = sm
        Kora._set_models!(gm, sm)
        _dhw_ref[] = nothing
        _init_n_ts_ref[] = 0
        _dhw_seed_ref[] = 0
        _dhw_override_ref[] = nothing
        _group_fraction_ref[] = (0.2f0, 0.2f0, 0.2f0, 0.2f0, 0.2f0)
        return Int32(0)
    catch e
        @_write_stderr("[bridge_aot] kf_load_models: ")
        @_write_exception(e)
        @_write_stderr("\n")
        return Int32(-1)
    end
end

Base.@ccallable function kf_set_deployment(
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

Base.@ccallable function kf_set_initial_cover(frac_ptr::Ptr{Float32}, n::Int32)::Int32
    n != Int32(5) && return Int32(-1)
    f1 = unsafe_load(frac_ptr, 1)
    f2 = unsafe_load(frac_ptr, 2)
    f3 = unsafe_load(frac_ptr, 3)
    f4 = unsafe_load(frac_ptr, 4)
    f5 = unsafe_load(frac_ptr, 5)
    _group_fraction_ref[] = (f1, f2, f3, f4, f5)
    return Int32(0)
end

Base.@ccallable function kf_set_dhw_trajectory(vals_ptr::Ptr{Float32}, n::Int32)::Int32
    if n <= Int32(0)
        _dhw_override_ref[] = nothing
        return Int32(0)
    end
    v = Vector{Float32}(undef, Int(n))
    GC.@preserve v unsafe_copyto!(pointer(v), vals_ptr, Int(n))
    _dhw_override_ref[] = v
    return Int32(0)
end

Base.@ccallable function kf_new_dhw_trajectory()::Int32
    _dhw_ref[] = nothing
    # A freshly requested seed trajectory supersedes any custom override.
    _dhw_override_ref[] = nothing
    return Int32(0)
end

Base.@ccallable function kf_run_reef(
    area_m2::Float32,
    init_cover_pct::Float32,
    n_runs::UInt32,
    dhw_seed::UInt32,
    dhw_out::Ptr{Float32},
    dhw_cap::Int32,
    covers_out::Ptr{Float32},
    covers_cap::Int32,
    lower_out::Ptr{Float32},
    median_out::Ptr{Float32},
    upper_out::Ptr{Float32},
    tol_lower_out::Ptr{Float32},
    tol_median_out::Ptr{Float32},
    tol_upper_out::Ptr{Float32},
    stats_cap::Int32,
    n_ts_out::Ptr{Int64},
    n_valid_out::Ptr{Int64}
)::Int32
    try
        dhw_cap <= Int32(0) && return Int32(-1)
        n_ts = Int(dhw_cap)

        gm = _growth_ref[]
        gm === nothing && return Int32(-1)
        sm = _survival_ref[]
        sm === nothing && return Int32(-1)

        # Generate DHW once per (n_ts, seed); reuse across run batches (and
        # across reef-area changes) so all runs see the same climate forcing
        # unless a new trajectory is explicitly requested via
        # kf_new_dhw_trajectory (which invalidates the cache; the seed for
        # the resulting regeneration is whatever this call passes in).
        _ov = _dhw_override_ref[]
        if _ov !== nothing && length(_ov) == n_ts
            @_write_stderr("[kf_run_reef] using custom DHW override\n")
            dhw_mat = reshape(copy(_ov), n_ts, 1)
        else
            if _dhw_ref[] === nothing || _init_n_ts_ref[] != n_ts || _dhw_seed_ref[] != dhw_seed
                @_write_stderr("[kf_run_reef] generate_example_dhw\n")
                _init_n_ts_ref[] = n_ts
                _dhw_seed_ref[] = dhw_seed
                _dhw_ref[] = Kora.generate_example_dhw(n_ts, 1; rng=Xoshiro(Int(dhw_seed)))
            end
            dhw_mat = _dhw_ref[]::Matrix{Float32}
        end
        unsafe_copyto!(dhw_out, pointer(dhw_mat[:, 1]), n_ts)

        @_write_stderr("[kf_run_reef] initialize_reef\n")
        reef = Kora.initialize_reef(;
            n_timesteps=n_ts, n_locs=1, area=Float64(area_m2), density=10,
            depths=Float64(_depth_ref[]),
            growth_models=gm, survival_models=sm
        )
        @_write_stderr("[kf_run_reef] initialize_coral_population!\n")
        Kora.initialize_coral_population!(reef)

        # Populate deployment schedule from globals set by kf_set_deployment.
        # reset!() does not clear deployment_times, so this persists per ensemble member.
        _apply_deployment!(reef, n_ts)

        dhw_tol = _deploy_dhw_tol_ref[]

        n_members = Int(n_runs)
        n_groups = Int(_N_GROUPS)
        ensemble_params = _build_ensemble_params(area_m2, init_cover_pct, n_members)

        @_write_stderr("[kf_run_reef] run_ensemble!\n")
        results = Kora.run_ensemble!(reef, dhw_mat, ensemble_params; deploy_dhw_tol=dhw_tol)
        @_write_stderr("[kf_run_reef] post-processing\n")

        valid_mask = [!any(isnan, results.cover[:, 1, r]) for r in 1:n_members]
        valid_indices = findall(valid_mask)
        n_valid = length(valid_indices)

        if covers_cap < n_ts * n_valid || stats_cap < n_ts * n_groups
            return Int32(-2)
        end

        for (col, r) in enumerate(valid_indices)
            for t in 1:n_ts
                unsafe_store!(covers_out, Float32(results.cover[t, 1, r]), (col-1)*n_ts + t)
            end
        end

        for g in 1:n_groups
            group_data = results.group_cover[:, 1, g, :]
            for t in 1:n_ts
                vals = filter(!isnan, collect(group_data[t, :]))
                lo, med, hi = if isempty(vals)
                    NaN32, NaN32, NaN32
                else
                    v = quantile(vals, [0.025, 0.5, 0.975])
                    Float32(v[1]), Float32(v[2]), Float32(v[3])
                end
                idx = (g-1)*n_ts + t
                unsafe_store!(lower_out, lo, idx)
                unsafe_store!(median_out, med, idx)
                unsafe_store!(upper_out, hi, idx)
            end
        end

        # Wild-population mean DHW tolerance per group per timestep. Same
        # percentile pass across members; dim 4 = [mean, std], take mean (index 1).
        for g in 1:n_groups
            tol_data = results.wild_dhw_tolerances[:, 1, g, 1, :]
            for t in 1:n_ts
                vals = filter(!isnan, collect(tol_data[t, :]))
                lo, med, hi = if isempty(vals)
                    NaN32, NaN32, NaN32
                else
                    v = quantile(vals, [0.025, 0.5, 0.975])
                    Float32(v[1]), Float32(v[2]), Float32(v[3])
                end
                idx = (g-1)*n_ts + t
                unsafe_store!(tol_lower_out, lo, idx)
                unsafe_store!(tol_median_out, med, idx)
                unsafe_store!(tol_upper_out, hi, idx)
            end
        end

        unsafe_store!(n_ts_out, Int64(n_ts))
        unsafe_store!(n_valid_out, Int64(n_valid))
        return Int32(n_ts)
    catch e
        @_write_stderr("[bridge_aot] kf_run_reef: ")
        @_write_exception(e)
        @_write_stderr("\n")
        return Int32(-1)
    end
end

end  # module KoraBridge
