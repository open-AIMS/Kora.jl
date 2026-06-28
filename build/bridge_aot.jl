# AOT entry points for kora_ui, compiled via:
#   juliac --project=julia_build --output-lib libkora_bridge --trim=safe bridge_aot.jl
#
# All @ccallable functions use C ABI: no Julia-managed types in signatures.
# Arrays are passed as (Ptr{T}, Int32 capacity); outputs are written in-place.
# The caller (Rust) must allocate output buffers before each call.
#
# State ordering: kf_init_reef MUST be called before kf_run_ensemble.

using Kora
using Statistics: quantile

# ── Module-level simulation state ─────────────────────────────────────────────

mutable struct _ReefState
    reef::Any
    dhw::Matrix{Float32}
    n_timesteps::Int
end

const _state = Ref{Union{Nothing,_ReefState}}(nothing)

const _N_GROUPS = Int32(5)
const _N_TIMESTEPS_DEFAULT = Int32(75)

# ── kf_init_reef ──────────────────────────────────────────────────────────────
#
# Initialise the reef and generate a synthetic DHW environment.
# Writes the DHW time series for location 1 into dhw_out (column-major Float32).
#
# Parameters:
#   area_m2        reef area in m²
#   init_cover_pct initial coral cover % (accepted; not yet wired into population)
#   dhw_out        caller-allocated Float32 buffer, capacity = dhw_cap
#   dhw_cap        capacity of dhw_out; must be >= n_timesteps
#
# Returns n_timesteps (75) on success, -1 on error or insufficient buffer.
Base.@ccallable function kf_init_reef(
    area_m2::Float32,
    ::Float32,                  # init_cover_pct — accepted for ABI stability; not yet wired
    dhw_out::Ptr{Float32},
    dhw_cap::Int32
)::Int32
    try
        n_ts = Int(_N_TIMESTEPS_DEFAULT)
        if dhw_cap < _N_TIMESTEPS_DEFAULT
            return Int32(-1)
        end

        reef = Kora.initialize_reef(;
            n_timesteps=n_ts,
            n_locs=1,
            area=Float64(area_m2),
            density=10,
            growth_models=Kora.growth_models,
            survival_models=Kora.survival_models
        )
        Kora.initialize_coral_population!(reef)
        dhw_mat::Matrix{Float32} = Kora.generate_example_dhw(n_ts, 1)

        _state[] = _ReefState(reef, dhw_mat, n_ts)

        dhw_loc1 = dhw_mat[:, 1]
        unsafe_copyto!(dhw_out, pointer(dhw_loc1), n_ts)

        return _N_TIMESTEPS_DEFAULT
    catch
        return Int32(-1)
    end
end

# ── kf_run_ensemble ───────────────────────────────────────────────────────────
#
# Run the ensemble using the current reef/env state and fill caller-allocated
# output buffers (all column-major Float32).
#
# Buffer layout (column-major):
#   covers_out  [n_timesteps × n_valid_runs]   total coral cover per run
#   lower_out   [n_timesteps × N_GROUPS]        2.5th-percentile group cover
#   median_out  [n_timesteps × N_GROUPS]        50th-percentile group cover
#   upper_out   [n_timesteps × N_GROUPS]        97.5th-percentile group cover
#
# Parameters:
#   deploy_volumes       Ptr to UInt32[n_vol]  (Phase 1: accepted, not applied)
#   n_vol                length of deploy_volumes
#   deploy_start_year    deployment start year (Phase 1: not applied)
#   deploy_cadence_years deployment cadence    (Phase 1: not applied)
#   n_runs               ensemble member count
#   covers_out / covers_cap   cover buffer and its capacity (Float32 elements)
#   lower_out, median_out, upper_out / stats_cap  stats buffers and capacity
#   n_ts_out             receives n_timesteps (Int64)
#   n_valid_out          receives n_valid_runs (Int64)
#
# Returns: 0 on success, -1 if uninitialised or exception, -2 if buffer too small.
Base.@ccallable function kf_run_ensemble(
    deploy_volumes::Ptr{UInt32},
    n_vol::Int32,
    deploy_start_year::UInt32,
    deploy_cadence_years::UInt32,
    n_runs::UInt32,
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
        state = _state[]
        state === nothing && return Int32(-1)

        n_ts = state.n_timesteps
        n_members = Int(n_runs)
        n_groups = Int(_N_GROUPS)

        ensemble_params = randn(6, n_members)
        ensemble_params[1, :] .= abs.(ensemble_params[1, :]) .+ 0.5
        for j in 1:n_members
            w = exp.(ensemble_params[2:6, j])
            ensemble_params[2:6, j] .= w ./ sum(w)
        end

        results = Kora.run_ensemble!(state.reef, state.dhw, ensemble_params)

        valid_mask = [!any(isnan, results.cover[:, 1, r]) for r in 1:n_members]
        valid_indices = findall(valid_mask)
        n_valid = length(valid_indices)

        if covers_cap < n_ts * n_valid || stats_cap < n_ts * n_groups
            return Int32(-2)
        end

        for (col, r) in enumerate(valid_indices)
            for t in 1:n_ts
                unsafe_store!(covers_out, Float32(results.cover[t, 1, r]),
                    (col - 1) * n_ts + t)
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
                idx = (g - 1) * n_ts + t
                unsafe_store!(lower_out, lo, idx)
                unsafe_store!(median_out, med, idx)
                unsafe_store!(upper_out, hi, idx)
            end
        end

        unsafe_store!(n_ts_out, Int64(n_ts))
        unsafe_store!(n_valid_out, Int64(n_valid))
        return Int32(0)
    catch
        return Int32(-1)
    end
end
