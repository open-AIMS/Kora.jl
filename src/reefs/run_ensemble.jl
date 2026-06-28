using TerminalLoggers, ProgressLogging

"""
    run_ensemble!(
        reef_state::ReefState,
        env_conditions::DimArray,
        ensemble_params::Matrix{Float64};
        rng::AbstractRNG=Random.GLOBAL_RNG
    )

Run an ensemble of simulations, one per column of `ensemble_params`, reusing
`reef_state` across members (reset between each run via `set_population!`).
Progress is reported to the terminal via `ProgressLogging`.

When `ensemble_params` has more than 16 rows, rows 17 through
`16 + n_groups(reef_state)` are interpreted as per-group location scalers passed
to `assign_scalers!`, and the final two rows as `recruits` and `self_seed`
parameters forwarded to `run_model!`.

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
    n_ensemble = size(ensemble_params, 2)
    n_ts = n_timesteps(reef_state)
    n_locs = n_locations(reef_state)
    n_grps = n_groups(reef_state)

    # Pre-allocate storage for ensemble results
    ensemble_cover = zeros(Float32, n_ts, n_locs, n_ensemble)
    ensemble_group_cover = zeros(Float32, n_ts, n_locs, n_grps, n_ensemble)
    ensemble_juvenile_cover = zeros(Float32, n_ts, n_locs, n_grps, n_ensemble)
    ensemble_wild_dhw_tolerances = zeros(Float32, n_ts, n_locs, n_grps, 2, n_ensemble)

    time_taken = @elapsed @progress for i in 1:n_ensemble
        params = ensemble_params[:, i]

        # Set up population with this parameter set
        set_population!(reef_state, params)

        # Apply scalers if present
        if length(params) > 16
            scaler_end = 17 + n_grps - 1
            loc_scalers = params[17:scaler_end]
            assign_scalers!(reef_state, loc_scalers)

            # Extract recruitment parameters
            recruitment_proportion = Float32(params[end - 1])
            self_seeding_proportion = Float32(params[end])

            # Run simulation
            run_model!(
                reef_state,
                env_conditions;
                recruits=recruitment_proportion,
                self_seed=self_seeding_proportion,
                rng=rng
            )
        else
            # Run with default recruitment parameters
            run_model!(reef_state, env_conditions; rng=rng)
        end

        mature_sizes = mature_size_thresholds()

        # Store results - optimized to avoid repeated function calls
        # and temporary allocations
        for ts in 1:n_ts
            for loc in 1:n_locs
                # Total cover at this timestep/location
                loc_cover = 0.0f0

                for grp in 1:n_grps
                    pop = coral_population(reef_state, ts, loc, grp)
                    grp_cover = sum(cover_cm_to_m2.(pop))

                    # Store group-level cover
                    ensemble_group_cover[ts, loc, grp, i] = grp_cover

                    ensemble_juvenile_cover[ts, loc, grp, i] = sum(
                        cover_cm_to_m2.(pop[pop .< mature_sizes[grp]])
                    )

                    # Store DHW tolerances (mean and stdev)
                    ensemble_wild_dhw_tolerances[ts, loc, grp, 1, i] = reef_state.wild_dhw_tolerances[
                        ts, loc, grp, 1
                    ]
                    ensemble_wild_dhw_tolerances[ts, loc, grp, 2, i] = reef_state.wild_dhw_tolerances[
                        ts, loc, grp, 2
                    ]

                    # Accumulate total cover
                    loc_cover += grp_cover
                end

                # Store total cover
                ensemble_cover[ts, loc, i] = loc_cover
            end
        end
    end

    @info "Ensemble completed in $(time_taken) seconds!"

    return (
        cover=ensemble_cover,
        group_cover=ensemble_group_cover,
        juvenile_cover=ensemble_juvenile_cover,
        wild_dhw_tolerances=ensemble_wild_dhw_tolerances,
        params=ensemble_params
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

    ensemble_cover = zeros(Float32, n_ts, n_locs, n_ensemble)
    ensemble_group_cover = zeros(Float32, n_ts, n_locs, n_grps, n_ensemble)
    ensemble_juvenile_cover = zeros(Float32, n_ts, n_locs, n_grps, n_ensemble)
    ensemble_wild_dhw_tolerances = zeros(Float32, n_ts, n_locs, n_grps, 2, n_ensemble)

    for i in 1:n_ensemble
        params = ensemble_params[:, i]
        set_population!(reef_state, params)

        if length(params) > 16
            scaler_end = 17 + n_grps - 1
            loc_scalers = params[17:scaler_end]
            assign_scalers!(reef_state, loc_scalers)
            recruitment_proportion = Float32(params[end - 1])
            self_seeding_proportion = Float32(params[end])
            run_model!(reef_state, dhw;
                recruits=recruitment_proportion,
                self_seed=self_seeding_proportion,
                rng=rng,
            )
        else
            run_model!(reef_state, dhw; rng=rng)
        end

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
                    ensemble_wild_dhw_tolerances[ts, loc, grp, 1, i] =
                        reef_state.wild_dhw_tolerances[ts, loc, grp, 1]
                    ensemble_wild_dhw_tolerances[ts, loc, grp, 2, i] =
                        reef_state.wild_dhw_tolerances[ts, loc, grp, 2]
                    loc_cover += grp_cover
                end
                ensemble_cover[ts, loc, i] = loc_cover
            end
        end
    end

    return (
        cover=ensemble_cover,
        group_cover=ensemble_group_cover,
        juvenile_cover=ensemble_juvenile_cover,
        wild_dhw_tolerances=ensemble_wild_dhw_tolerances,
        params=ensemble_params
    )
end

"""
    run_ensemble!(
        reef_state::ReefState,
        env_conditions::DimArray,
        ensemble_params::Matrix{Float64},
        ::Val{:extended};
        rng::AbstractRNG=Random.GLOBAL_RNG
    )

Extended overload that always applies location scalers and recruitment parameters.
Intended for parameter matrices with more than 16 rows.

Expected row layout per column (n = `n_groups(reef_state)`):

| Rows     | Content                                        |
|----------|------------------------------------------------|
| 1–16     | Population parameters consumed by `set_population!` |
| 17–16+n  | Per-group growth scalers passed to `assign_scalers!` |
| end-1    | `recruits` proportion forwarded to `run_model!`  |
| end      | `self_seed` proportion forwarded to `run_model!` |

For a 5-group model this requires 16 + 5 + 2 = **23 rows**.

# See Also
[`run_ensemble!`](@ref), [`run_model!`](@ref), [`set_population!`](@ref), [`assign_scalers!`](@ref)
"""
function run_ensemble!(
    reef_state::ReefState,
    env_conditions::DimArray,
    ensemble_params::Matrix{Float64},
    ::Val{:extended};
    rng::AbstractRNG=Random.GLOBAL_RNG
)
    n_ensemble = size(ensemble_params, 2)
    n_ts = n_timesteps(reef_state)
    n_locs = n_locations(reef_state)
    n_grps = n_groups(reef_state)

    n_params = size(ensemble_params, 1)
    expected = 16 + n_grps + 2
    if n_params != expected
        throw(ArgumentError(
            "ensemble_params has $n_params rows; expected $expected " *
            "(16 population + $n_grps scalers + 2 recruitment)"
        ))
    end

    scaler_end = 16 + n_grps

    ensemble_cover = zeros(Float32, n_ts, n_locs, n_ensemble)
    ensemble_group_cover = zeros(Float32, n_ts, n_locs, n_grps, n_ensemble)
    ensemble_juvenile_cover = zeros(Float32, n_ts, n_locs, n_grps, n_ensemble)
    ensemble_wild_dhw_tolerances = zeros(Float32, n_ts, n_locs, n_grps, 2, n_ensemble)

    time_taken = @elapsed @progress for i in 1:n_ensemble
        params = ensemble_params[:, i]

        set_population!(reef_state, params)

        loc_scalers = params[17:scaler_end]
        assign_scalers!(reef_state, loc_scalers)

        recruitment_proportion = Float32(params[end - 1])
        self_seeding_proportion = Float32(params[end])

        run_model!(
            reef_state,
            env_conditions;
            recruits=recruitment_proportion,
            self_seed=self_seeding_proportion,
            rng=rng
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

    @info "Ensemble completed in $(time_taken) seconds!"

    return (
        cover=ensemble_cover,
        group_cover=ensemble_group_cover,
        juvenile_cover=ensemble_juvenile_cover,
        wild_dhw_tolerances=ensemble_wild_dhw_tolerances,
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
