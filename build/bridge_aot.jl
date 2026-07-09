# bridge_aot.jl — Julia AOT bridge for Kora reef simulation
#
# C API signatures (for Rust FFI declarations):
#   int32_t kf_load_models(const uint8_t* growth_path, const uint8_t* surv_path);
#   int32_t kf_set_deployment(uint32_t vol0, uint32_t vol1, uint32_t vol2,
#                             uint32_t vol3, uint32_t vol4,
#                             uint32_t start_year, uint32_t cadence_years,
#                             float depth_m, float deploy_dhw_tolerance);
#   int32_t kf_run_reef(float area_m2, float init_cover_pct, uint32_t n_runs,
#                       float* dhw_out, int32_t dhw_cap,
#                       float* covers_out, int32_t covers_cap,
#                       float* lower_out, float* median_out, float* upper_out,
#                       int32_t stats_cap,
#                       int64_t* n_ts_out, int64_t* n_valid_out);

module KoraBridge

using Kora
using Statistics: quantile

macro _write_stderr(msg)
    n = ncodeunits(msg)
    :(ccall(:_write, Cint, (Cint, Ptr{UInt8}, Cuint), Int32(2), $msg, Cuint($n)))
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
            ccall(:_write, Cint, (Cint, Ptr{UInt8}, Cuint), Int32(2), _tcstr, Cuint(_tlen))
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
        ccall(:_write, Cint, (Cint, Ptr{UInt8}, Cuint), Int32(2), " on ", Cuint(4))
        ccall(:_write, Cint, (Cint, Ptr{UInt8}, Cuint), Int32(2), fstr, Cuint(flen))
    elseif exc isa CompositeException
        excs = (exc::CompositeException).exceptions
        if !isempty(excs)
            tfe = excs[1]
            if tfe isa TaskFailedException
                inner = (tfe::TaskFailedException).task.result
                icstr = ccall(:jl_typeof_str, Ptr{UInt8}, (Any,), inner)
                ilen = ccall(:strlen, Csize_t, (Ptr{UInt8},), icstr)
                ccall(
                    :_write,
                    Cint,
                    (Cint, Ptr{UInt8}, Cuint),
                    Int32(2),
                    " [inner: ",
                    Cuint(9)
                )
                ccall(
                    :_write, Cint, (Cint, Ptr{UInt8}, Cuint), Int32(2), icstr, Cuint(ilen)
                )
                if inner isa MethodError
                    fn = (inner::MethodError).f
                    fstr = ccall(:jl_typeof_str, Ptr{UInt8}, (Any,), fn)
                    flen = ccall(:strlen, Csize_t, (Ptr{UInt8},), fstr)
                    ccall(
                        :_write, Cint, (Cint, Ptr{UInt8}, Cuint), Int32(2), " on ", Cuint(4)
                    )
                    ccall(
                        :_write,
                        Cint,
                        (Cint, Ptr{UInt8}, Cuint),
                        Int32(2),
                        fstr,
                        Cuint(flen)
                    )
                end
                ccall(:_write, Cint, (Cint, Ptr{UInt8}, Cuint), Int32(2), "]", Cuint(1))
            end
        end
    end
    return nothing
end

const _N_TIMESTEPS_DEFAULT = Int32(50)
const _N_GROUPS = Int32(5)

const _growth_ref = Ref{Union{Nothing,Kora.PolyGrowthModel{Float32}}}(nothing)
const _survival_ref = Ref{Union{Nothing,Kora.PolySurvivalModel{Float32}}}(nothing)

# Cache DHW so every run batch uses the same climate forcing.
# Regenerated when reef area changes or models are reloaded.
const _dhw_ref = Ref{Union{Nothing,Matrix{Float32}}}(nothing)
const _init_area_ref = Ref{Float32}(0.0f0)

# Deployment schedule — set via kf_set_deployment before kf_run_reef.
# NTuple{5,UInt32}: corals/year for each of the 5 functional groups.
const _deploy_vols_ref = Ref{NTuple{5,UInt32}}((
    UInt32(0), UInt32(0), UInt32(0), UInt32(0), UInt32(0)
))
const _deploy_start_ref = Ref{UInt32}(UInt32(1))
const _deploy_cadence_ref = Ref{UInt32}(UInt32(1))
const _depth_ref = Ref{Float32}(9.0f0)
const _deploy_dhw_tol_ref = Ref{Float32}(0.0f0)

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
    params[2:6, :] .= 0.2
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
        _init_area_ref[] = 0.0f0
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

Base.@ccallable function kf_run_reef(
    area_m2::Float32,
    init_cover_pct::Float32,
    n_runs::UInt32,
    dhw_out::Ptr{Float32},
    dhw_cap::Int32,
    covers_out::Ptr{Float32},
    covers_cap::Int32,
    lower_out::Ptr{Float32},
    median_out::Ptr{Float32},
    upper_out::Ptr{Float32},
    stats_cap::Int32,
    n_ts_out::Ptr{Int64},
    n_valid_out::Ptr{Int64}
)::Int32
    try
        n_ts = Int(_N_TIMESTEPS_DEFAULT)
        dhw_cap < _N_TIMESTEPS_DEFAULT && return Int32(-1)

        gm = _growth_ref[]
        gm === nothing && return Int32(-1)
        sm = _survival_ref[]
        sm === nothing && return Int32(-1)

        # Generate DHW once per area; reuse across run batches so all runs see the same
        # climate forcing.  Cleared by kf_load_models when models are reloaded.
        if _dhw_ref[] === nothing || _init_area_ref[] != area_m2
            @_write_stderr("[kf_run_reef] generate_example_dhw\n")
            _init_area_ref[] = area_m2
            _dhw_ref[] = Kora.generate_example_dhw(n_ts, 1)
        end
        dhw_mat = _dhw_ref[]::Matrix{Float32}
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
        if dhw_tol > 0.0f0
            dhw_tol_sigma = 0.5f0   # fixed narrow spread; adjust if calibration data is available
            vols = _deploy_vols_ref[]
            start = Int(_deploy_start_ref[])
            cadence = Int(_deploy_cadence_ref[])
            for ts in start:cadence:n_ts, grp in 1:5
                if vols[grp] > 0
                    reef.deployed_dhw_tolerances[ts, 1, grp, 1] = dhw_tol
                    reef.deployed_dhw_tolerances[ts, 1, grp, 2] = dhw_tol_sigma
                end
            end
        end

        n_members = Int(n_runs)
        n_groups = Int(_N_GROUPS)
        ensemble_params = _build_ensemble_params(area_m2, init_cover_pct, n_members)

        @_write_stderr("[kf_run_reef] run_ensemble!\n")
        results = Kora.run_ensemble!(reef, dhw_mat, ensemble_params)
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

        unsafe_store!(n_ts_out, Int64(n_ts))
        unsafe_store!(n_valid_out, Int64(n_valid))
        return _N_TIMESTEPS_DEFAULT
    catch e
        @_write_stderr("[bridge_aot] kf_run_reef: ")
        @_write_exception(e)
        @_write_stderr("\n")
        return Int32(-1)
    end
end

end  # module KoraBridge
