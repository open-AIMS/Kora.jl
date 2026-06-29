function _alloc_ensemble_results(
    n_ts::Int, n_locs::Int, n_grps::Int, n_ensemble::Int
)
    return (
        zeros(Float32, n_ts, n_locs, n_ensemble),
        zeros(Float32, n_ts, n_locs, n_grps, n_ensemble),
        zeros(Float32, n_ts, n_locs, n_grps, n_ensemble),
        zeros(Float32, n_ts, n_locs, n_grps, 2, n_ensemble)
    )
end

function _collect_member!(
    ensemble_cover::Array{Float32,3},
    ensemble_group_cover::Array{Float32,4},
    ensemble_juvenile_cover::Array{Float32,4},
    ensemble_wild_dhw_tolerances::Array{Float32,5},
    reef_state::ReefState,
    i::Int,
    n_ts::Int,
    n_locs::Int,
    n_grps::Int
)
    mature_sizes = mature_size_thresholds()
    for ts in 1:n_ts
        for loc in 1:n_locs
            loc_cover = 0.0f0
            for grp in 1:n_grps
                pop = coral_population(reef_state, ts, loc, grp)
                grp_cover = sum(cover_cm_to_m2.(pop))
                ensemble_group_cover[ts, loc, grp, i] = grp_cover
                ensemble_juvenile_cover[ts, loc, grp, i] = sum(
                    cover_cm_to_m2.(pop[pop .< mature_sizes[grp]])
                )
                ensemble_wild_dhw_tolerances[ts, loc, grp, 1, i] = reef_state.wild_dhw_tolerances[
                    ts, loc, grp, 1
                ]
                ensemble_wild_dhw_tolerances[ts, loc, grp, 2, i] = reef_state.wild_dhw_tolerances[
                    ts, loc, grp, 2
                ]
                loc_cover += grp_cover
            end
            ensemble_cover[ts, loc, i] = loc_cover
        end
    end
end

"""
    run_ensemble!(
        reef_state::ReefState,
        env_conditions::DimArray,
        ensemble_params::Matrix{Float64};
        rng::AbstractRNG=Random.GLOBAL_RNG
    )

Run an ensemble of simulations, one per column of `ensemble_params`, reusing
`reef_state` across members (reset between each run via `set_population!`).

When `ensemble_params` has more than 16 rows, rows 17 through
`16 + n_groups(reef_state)` are interpreted as per-group location scalers passed
to `assign_scalers!`, and the final two rows as `recruits` and `self_seed`
parameters forwarded to `run_model!`.

`env_conditions` must be a 3D `DimArray` with axes
`(Dim{:timestep}, Dim{:location}, Dim{:variable})` containing at minimum a
`:dhw` variable slice.

# Arguments
- `reef_state` : `ReefState` used as the simulation template. Mutated during
  each member run; contents after the call reflect only the last ensemble member.
- `env_conditions` : Environmental forcing data (DimArray) shared across all members.
- `ensemble_params` : Parameter matrix of shape `(n_params, n_members)`. Each
  column defines one ensemble member. Rows 1-16 are population parameters
  consumed by `set_population!`.
- `rng` : Random number generator (default: `Random.GLOBAL_RNG`).

# Returns
`NamedTuple` with the fields listed below.

- `cover::Array{Float32,3}` : Total coral cover in m^2 with shape
  `(n_timesteps, n_locations, n_members)`.
- `group_cover::Array{Float32,4}` : Per-group cover in m^2 with shape
  `(n_timesteps, n_locations, n_groups, n_members)`.
- `juvenile_cover::Array{Float32,4}` : Sub-mature coral cover in m^2 with shape
  `(n_timesteps, n_locations, n_groups, n_members)`.
- `wild_dhw_tolerances::Array{Float32,5}` : DHW tolerance statistics with shape
  `(n_timesteps, n_locations, n_groups, 2, n_members)`. The third inner
  dimension holds mean (index 1) and standard deviation (index 2).
- `params::Matrix{Float64}` : The input `ensemble_params` unchanged.

# See Also
[`run_model!`](@ref), [`set_population!`](@ref), [`assign_scalers!`](@ref)
"""
function run_ensemble!(
    reef_state::ReefState,
    env_conditions::DimArray,
    ensemble_params::Matrix{Float64};
    rng::AbstractRNG=Random.GLOBAL_RNG
)
    return run_ensemble!(
        reef_state,
        Matrix{Float32}(env_conditions[:, :, At(:dhw)].data),
        ensemble_params;
        rng=rng
    )
end

function run_ensemble!(
    reef_state::ReefState,
    dhw::Matrix{Float32},
    ensemble_params::Matrix{Float64};
    rng::AbstractRNG=Random.GLOBAL_RNG
)
    n_ensemble = size(ensemble_params, 2)
    n_ts = n_timesteps(reef_state)
    n_locs = n_locations(reef_state)
    n_grps = n_groups(reef_state)

    ec, egc, ejc, ewdt = _alloc_ensemble_results(n_ts, n_locs, n_grps, n_ensemble)

    for i in 1:n_ensemble
        params = ensemble_params[:, i]
        set_population!(reef_state, params)
        if length(params) > 16
            expected = 16 + n_grps + 2
            length(params) == expected || throw(ArgumentError(
                "ensemble_params has $(length(params)) rows; expected $expected " *
                "(16 population + $n_grps scalers + 2 recruitment)"
            ))
            assign_scalers!(reef_state, params[17:(16 + n_grps)])
            run_model!(reef_state, dhw;
                recruits=Float32(params[end - 1]),
                self_seed=Float32(params[end]),
                rng=rng
            )
        else
            run_model!(reef_state, dhw; rng=rng)
        end
        _collect_member!(ec, egc, ejc, ewdt, reef_state, i, n_ts, n_locs, n_grps)
    end

    return (
        cover=ec,
        group_cover=egc,
        juvenile_cover=ejc,
        wild_dhw_tolerances=ewdt,
        params=ensemble_params
    )
end

"""
    summarize_ensemble(ensemble_results, area; quantiles=[0.025, 0.5, 0.975])

Calculate summary statistics across ensemble members.

# Returns
NamedTuple with mean, median, and quantiles of ensemble predictions
"""
function summarize_ensemble(ensemble_results, area; quantiles=[0.025, 0.5, 0.975])
    cover = ensemble_results.cover
    n_ts, n_locs, n_ensemble = size(cover)

    # Convert to percentage
    cover_pct = (cover ./ area) .* 100.0

    # Calculate statistics across ensemble dimension
    mean_cover = mean(cover_pct; dims=3)[:, :, 1]
    median_cover = median(cover_pct; dims=3)[:, :, 1]

    # Calculate quantiles
    quantile_cover = zeros(Float32, n_ts, n_locs, length(quantiles))
    for ts in 1:n_ts, loc in 1:n_locs
        quantile_cover[ts, loc, :] = quantile(cover_pct[ts, loc, :], quantiles)
    end

    return (
        mean=mean_cover,
        median=median_cover,
        quantiles=quantile_cover,
        quantile_levels=quantiles
    )
end
