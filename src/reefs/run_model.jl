"""
    run_model!(
        reef_state::ReefState,
        env_conditions::YAXArray;
        recruits=0.06f0,
        self_seed=0.3f0,
        rng::AbstractRNG=Random.GLOBAL_RNG
    )::Nothing

Advance the reef simulation forward in time, writing results into `reef_state`
in-place. The state is reset to its timestep-1 population before the run begins,
so calling `run_model!` a second time on the same object produces a fresh result.

`env_conditions` must be a 3D `YAXArray` with axes
`(Dim{:timestep}, Dim{:location}, Dim{:variable})` containing at minimum a
`:dhw` variable slice. This is the format returned by both
`generate_example_environment` and `generate_environment`.

# Arguments
- `reef_state` : Pre-initialised `ReefState`. Modified in-place; the caller's
  object contains the full time series after this call returns.
- `env_conditions` : Environmental forcing data. Must cover the same number of
  timesteps and locations as `reef_state`.
- `recruits` : Fraction of local larval production that successfully recruits to
  the reef each timestep (default: `0.06`).
- `self_seed` : Fraction of recruitment attributed to self-seeding from the
  local population (default: `0.3`).
- `rng` : Random number generator. Pass a seeded `Xoshiro` or similar for
  reproducible runs (default: `Random.GLOBAL_RNG`).

# Returns
`Nothing`

# See Also
[`initialize_reef`](@ref), [`initialize_coral_population!`](@ref),
[`run_model`](@ref), [`coral_cover`](@ref)
"""
function run_model!(
    reef_state::ReefState,
    env_conditions::YAXArray;
    recruits=0.06f0,
    self_seed=0.3f0,
    rng::AbstractRNG=Random.GLOBAL_RNG
)::Nothing
    reset!(reef_state)

    timesteps::UnitRange{Int64} = 1:n_timesteps(reef_state)
    n_locs::Int64 = n_locations(reef_state)
    n_grps::Int64 = n_groups(reef_state)

    # Assumed proportion of larvae contributing to coral recruitment
    recruitment_proportion = recruits
    self_seeding_proportion = self_seed

    # Deployment per group
    # Assuming a 200m study area, this comes to ~7 deployments per m²
    # which is relatively high deployment density
    # (n * n_grps) / area
    # (280 * 5) / 200 = 7
    n_deploy = 0
    recruit_μ = 1.5
    recruit_σ = 0.2

    # Mock recruitment distribution
    # TODO: Compare/contrast with other approaches
    recruit_dist = truncated(Normal(recruit_μ, recruit_σ), 0.5, 2.5)

    # Convenience var
    carrying_cap = reef_state.carrying_capacity
    total_possible_colonies = carrying_cap .* reef_state.density

    # Group modifiers that adjusts when growth slows down as it approaches a
    # proportionate limit relative to available space
    inflection_points = growth_inflection_point()
    maturity_thresholds = mature_size_thresholds()

    recruits = fill(Float32[], n_locs, n_grps)
    depth_coeffs = depth_coefficient.(reef_state.depths)
    pop_buffer = reef_state._pop_buffer

    # Clear any existing results
    reset!(reef_state)

    # Apply initial bleaching mortality (for ts = 1)
    # TODO: Abstract into separate callable method
    dhws = env_conditions[1, :, At(:dhw)].data
    for loc in 1:n_locs, grp in 1:n_grps
        pop_buffer .= 0.0f0  # Reset buffer

        fill_population_buffer!(
            reef_state, 1, loc, grp, recruits[loc, grp], pop_buffer
        )

        # Diameters for entire population including recruits
        with_recruits = pop_buffer[pop_buffer .> 0.0f0]

        # Background mortality
        # TODO: apply location specific survival scaler
        apply_survival!(reef_state, grp, with_recruits, rng)

        # Bleaching mortality
        tols = @view(reef_state.wild_dhw_tolerances.data[1, loc, grp, :])
        new_mean, new_std, area_lost = bleaching_mortality!(
            with_recruits,
            dhws[loc],
            depth_coeffs[loc],
            tols,
            grp
        )
        tols[1] = new_mean
        tols[2] = new_std

        # Cyclone mortality
        # cyclone_probs = cyclone_mortality_prob.(p_sample, [cyclone_cat])
        # total_probs = cyclone_probs  # TODO: Add other probabilities
        # survivors_mask = rand(length(p_sample)) .> total_probs

        update_pop_cache!(reef_state, with_recruits, loc)

        next_pop = @view(reef_state._pop_cache[loc, :])
        update_wild_sample!(reef_state, 1, loc, grp, next_pop[next_pop .> 0.0])
    end

    recruits = fill(Float32[], n_locs, n_grps)  # Create cache

    # TODO: Refactor into `run_timestep()`
    for ts in timesteps[2:end]
        prev_ts = ts - 1

        # Used to determine available space, which affects settlement rate, etc.
        total_covers::Vector{Float32} = coral_cover(reef_state, prev_ts)

        # Reset cache
        for i in eachindex(recruits)
            empty!(recruits[i])
        end

        for loc in 1:n_locs, grp in 1:n_grps
            # Some recruitment happens
            # Calculate across all locations and groups for time `ts`...

            # Scale recruitment by cover proportion
            prod = larval_production(reef_state, maturity_thresholds, prev_ts, loc, grp)
            prod = prod * self_seeding_proportion * recruitment_proportion

            # Only allow recruitment if density threshold has not been reached
            # Note this is natural recruitment. Deployments are handled separately.
            all_pop = total_population(reef_state, prev_ts, loc)
            avail_d_for_recruitment = all_pop < total_possible_colonies[loc]
            if avail_d_for_recruitment
                available_space = max(carrying_cap[loc] - total_covers[loc], 0.0)
                prod = (prod / reef_state.carrying_capacity[loc]) * available_space
                # n% of produced larvae arrive and a proportion (based on available area)
                # of these can settle
                settlement_prop = min.(
                    (available_space * 50),
                    prod
                )
                n_loc_recruits = floor(Int64, settlement_prop)
            else
                n_loc_recruits = 0
            end

            if n_loc_recruits > 0
                # Only selection changes tolerance
                recruits[loc, grp] = Float32.(rand(rng, recruit_dist, n_loc_recruits))
                update_coral_tolerances!(reef_state, ts, loc, grp, n_loc_recruits)
            else
                # Copy previous tolerance when no recruits (no Breeder's equation applied)
                reef_state.wild_dhw_tolerances.data[ts, loc, grp, 1] = reef_state.wild_dhw_tolerances.data[
                    prev_ts, loc, grp, 1
                ]
                reef_state.deployed_dhw_tolerances.data[ts, loc, grp, 1] = reef_state.deployed_dhw_tolerances.data[
                    prev_ts, loc, grp, 1
                ]
            end

            # Deploy scheduled corals (outplanted corals are mature and can reproduce immediately)
            if reef_state.deployment_times[ts, loc, grp] > 0
                n_deploy = Int64(reef_state.deployment_times[ts, loc, grp])
                deploy_corals!(reef_state, ts, loc, n_deploy, grp; rng=rng)
            end
        end

        dhws = env_conditions[ts, :, At(:dhw)].data
        # cyclone_cats = env_conditions[ts, :, At(:cyclone_category)].data

        # Mortality occurs
        for loc in 1:n_locs, grp in 1:n_grps
            pop_buffer .= 0.0f0  # Reset buffer

            fill_population_buffer!(
                reef_state, prev_ts, loc, grp, recruits[loc, grp], pop_buffer
            )

            # Diameters for entire population including recruits
            with_recruits = pop_buffer[pop_buffer .> 0.0f0]

            # Background mortality
            # TODO: apply location specific survival scaler
            apply_survival!(reef_state, grp, with_recruits, rng)

            # Bleaching mortality
            tols = @view(reef_state.wild_dhw_tolerances.data[ts, loc, grp, :])
            new_mean, new_std, area_lost = bleaching_mortality!(
                with_recruits,
                dhws[loc],
                depth_coeffs[loc],
                tols,
                grp
            )
            tols[1] = new_mean
            tols[2] = new_std

            # Cyclone mortality
            # cyclone_probs = cyclone_mortality_prob.(p_sample, [cyclone_cat])
            # total_probs = cyclone_probs  # TODO: Add other probabilities
            # survivors_mask = rand(length(p_sample)) .> total_probs

            update_pop_cache!(reef_state, with_recruits, loc)

            next_pop = @view(reef_state._pop_cache[loc, :])
            update_wild_sample!(reef_state, ts, loc, grp, next_pop[next_pop .> 0.0])
        end

        # Take snapshot of cover prior to growth
        total_covers = coral_cover(reef_state, ts - 1)
        cover_proportions = total_covers ./ reef_state.carrying_capacity

        # Survivers grow...
        # If any location is near capacity, refresh cover_proportions after each group
        # so that fast-growing groups (processed first, 1→5: Acropora → massives) consume
        # space before slower groups — preventing within-timestep overshoot.
        if any(>(0.85f0), cover_proportions)
            for grp in 1:n_grps
                apply_growth!(
                    reef_state,
                    grp,
                    inflection_points[grp],
                    reef_state.wild_population[ts, :, grp],
                    cover_proportions
                )
                total_covers = coral_cover(reef_state, ts)
                cover_proportions = total_covers ./ reef_state.carrying_capacity

                apply_growth!(
                    reef_state,
                    grp,
                    inflection_points[grp],
                    reef_state.deployed_population[ts, :, grp],
                    cover_proportions
                )
                total_covers = coral_cover(reef_state, ts)
                cover_proportions = total_covers ./ reef_state.carrying_capacity
            end
        else
            for grp in 1:n_grps
                apply_growth!(
                    reef_state,
                    grp,
                    inflection_points[grp],
                    reef_state.wild_population[ts, :, grp],
                    cover_proportions
                )
                apply_growth!(
                    reef_state,
                    grp,
                    inflection_points[grp],
                    reef_state.deployed_population[ts, :, grp],
                    cover_proportions
                )
            end
            total_covers = coral_cover(reef_state, ts)
        end

        # The sigmoid space_constraint can still allow small within-step overshoot.
        # Proportionally scale down diameters at any location that exceeded K,
        # preserving the relative size structure (cover ∝ d², so scale d by √(K/C)).
        # A small conservative margin (8 ulps) on the scale factor absorbs Float32
        # summation rounding that would otherwise leave recomputed cover 1 ulp above K.
        total_covers = coral_cover(reef_state, ts)
        for loc in 1:n_locs
            if total_covers[loc] > reef_state.carrying_capacity[loc]
                scale = sqrt(reef_state.carrying_capacity[loc] / total_covers[loc])
                scale *= 1.0f0 - 8.0f0 * eps(Float32)
                for grp in 1:n_grps
                    reef_state.wild_population[ts, loc, grp] .*= scale
                    reef_state.deployed_population[ts, loc, grp] .*= scale
                end
            end
        end
    end

    return nothing
end

"""
    run_model(;
        n_ts=75,
        n_locs=100,
        with_dhw=true,
        area=100.0,
        pop_density=15.0,
        growth_models=growth_models,
        survival_models=survival_models
    )::Tuple{ReefState, YAXArray}

Convenience wrapper that allocates a reef, seeds its population, generates
synthetic environmental conditions, and returns the completed simulation results.

Internally calls `initialize_reef`, `initialize_coral_population!`,
`generate_example_environment`, and `run_model!` in sequence using the supplied
keyword arguments.

# Arguments
- `n_ts` : Number of annual time steps (default: `75`).
- `n_locs` : Number of reef locations (default: `100`).
- `with_dhw` : Whether to generate DHW thermal forcing. Pass `false` to run
  with zero thermal stress (default: `true`).
- `area` : Reef area in m^2 used for carrying capacity (default: `100.0`).
- `pop_density` : Initial colony density in colonies per m^2 used to size the
  starting population (default: `15.0`).
- `growth_models` : Fitted growth model collection (default: package-level
  offshore-north models).
- `survival_models` : Fitted survival model collection (default: package-level
  offshore-north models).

# Returns
`Tuple{ReefState, YAXArray}` : The completed reef state containing the full
time series and the environmental conditions used for the run.

# See Also
[`run_model!`](@ref), [`initialize_reef`](@ref),
[`generate_example_environment`](@ref)
"""
function run_model(;
    n_ts=75, n_locs=100, with_dhw=true, area=100.0, pop_density=15.0,
    growth_models=growth_models, survival_models=survival_models
)
    sample_size = ceil(Int64, area * pop_density) * 2
    reef_state = initialize_reef(;
        n_timesteps=n_ts,
        n_locs=n_locs,
        area=area,
        density=pop_density,
        sample_size=sample_size,
        growth_models=growth_models,
        survival_models=survival_models
    )
    initialize_coral_population!(reef_state)
    env_template = generate_example_environment(n_ts, n_locs; with_dhw=with_dhw)

    run_model!(reef_state, env_template)

    return reef_state, env_template
end
