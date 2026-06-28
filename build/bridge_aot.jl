# bridge_aot.jl — Julia AOT bridge for Kora reef simulation
#
# C API signatures (for Rust FFI declarations):
#   int32_t kf_load_models(const uint8_t* growth_path, const uint8_t* surv_path);
#   int32_t kf_run_reef(float area_m2, float _reserved, uint32_t n_runs,
#                       float* dhw_out, int32_t dhw_cap,
#                       float* covers_out, int32_t covers_cap,
#                       float* lower_out, float* median_out, float* upper_out,
#                       int32_t stats_cap,
#                       int64_t* n_ts_out, int64_t* n_valid_out);

module KoraBridge

using Kora
using Statistics: quantile

const _N_TIMESTEPS_DEFAULT = Int32(50)
const _N_GROUPS = Int32(6)

const _growth_ref   = Ref{Union{Nothing, Kora.PolyGrowthModel}}(nothing)
const _survival_ref = Ref{Union{Nothing, Kora.PolySurvivalModel}}(nothing)

function _build_ensemble_params(n_members::Int)::Matrix{Float64}
    params = randn(6, n_members)
    params[1, :] .= abs.(params[1, :]) .+ 0.5
    for j in 1:n_members
        w = exp.(params[2:6, j])
        params[2:6, j] .= w ./ sum(w)
    end
    return params
end

Base.@ccallable function kf_load_models(
    growth_path::Ptr{UInt8},
    surv_path::Ptr{UInt8}
)::Int32
    try
        gp = unsafe_string(growth_path)
        sp = unsafe_string(surv_path)
        gm = Kora.load_models(gp)::Kora.PolyGrowthModel
        sm = Kora.load_models(sp)::Kora.PolySurvivalModel
        _growth_ref[]   = gm
        _survival_ref[] = sm
        Kora._set_models!(gm, sm)
        return Int32(0)
    catch
        return Int32(-1)
    end
end

Base.@ccallable function kf_run_reef(
    area_m2::Float32,
    ::Float32,
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

        reef = Kora.initialize_reef(;
            n_timesteps=n_ts, n_locs=1, area=Float64(area_m2), density=10,
            growth_models=gm, survival_models=sm
        )
        Kora.initialize_coral_population!(reef)

        dhw_mat = Kora.generate_example_dhw(n_ts, 1)
        unsafe_copyto!(dhw_out, pointer(dhw_mat[:, 1]), n_ts)

        n_members = Int(n_runs)
        n_groups  = Int(_N_GROUPS)
        ensemble_params = _build_ensemble_params(n_members)

        results = Kora.run_ensemble!(reef, dhw_mat, ensemble_params)

        valid_mask    = [!any(isnan, results.cover[:, 1, r]) for r in 1:n_members]
        valid_indices = findall(valid_mask)
        n_valid       = length(valid_indices)

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
                unsafe_store!(lower_out,  lo,  idx)
                unsafe_store!(median_out, med, idx)
                unsafe_store!(upper_out,  hi,  idx)
            end
        end

        unsafe_store!(n_ts_out,    Int64(n_ts))
        unsafe_store!(n_valid_out, Int64(n_valid))
        return _N_TIMESTEPS_DEFAULT
    catch
        return Int32(-1)
    end
end

end  # module KoraBridge
