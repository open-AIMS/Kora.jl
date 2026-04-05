using TerminalLoggers, ProgressLogging

"""
    run_ensemble!(
        reef_state::ReefState,
        env_conditions::YAXArray,
        ensemble_params::Matrix{Float64};
        rng::AbstractRNG=Random.GLOBAL_RNG
    )

Run ensemble simulations with multiple parameter sets.

# Arguments
- `reef_state` : ReefState (will be reused for each ensemble member)
- `env_conditions` : Environmental conditions
- `ensemble_params` : Matrix where each column is a parameter set
- `rng` : Random number generator state

# Returns
- `ensemble_results` : NamedTuple with results from all ensemble members
"""
function run_ensemble!(
    reef_state::ReefState,
    env_conditions::YAXArray,
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

    @info "Running ensemble of $(n_ensemble) simulations..."

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
            run_example!(
                reef_state,
                env_conditions;
                recruits=recruitment_proportion,
                self_seed=self_seeding_proportion,
                rng=rng
            )
        else
            # Run with default recruitment parameters
            run_example!(reef_state, env_conditions; rng=rng)
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
                        ts, loc, grp, At(:mean)
                    ]
                    ensemble_wild_dhw_tolerances[ts, loc, grp, 2, i] = reef_state.wild_dhw_tolerances[
                        ts, loc, grp, At(:stdev)
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
