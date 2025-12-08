"""
    calculate_inheritance_proportions(reef_state::ReefState, ts::Int64, loc::Int64, grp::Int64)::Tuple{Float32, Float32}

Calculate mixing proportions between wild and deployed populations based on their relative cover.

Returns (wild_proportion, deployed_proportion)
"""
function calculate_inheritance_proportions(
    reef_state::ReefState, ts::Int64, loc::Int64, grp::Int64
)::Tuple{Float32,Float32}
    # Calculate cover for each population
    wild_pop = wild_population(reef_state, ts, loc, grp)
    deployed_pop = deployed_population(reef_state, ts, loc, grp)

    wild_cover = sum(cover_cm_to_m2.(wild_pop))
    deployed_cover = sum(cover_cm_to_m2.(deployed_pop))
    total_cover = wild_cover + deployed_cover

    if total_cover == 0
        return (1.0f0, 0.0f0)  # Default to wild if no cover
    end

    return (wild_cover / total_cover, deployed_cover / total_cover)
end

"""
    update_coral_tolerances!(
        reef_state::ReefState,
        ts::Int64,
        loc::Int64,
        grp::Int64;
        h²::Float32=0.3f0
    )::Nothing

Update thermal tolerances for new recruits based on mixing between wild and deployed populations.

# Arguments
- `reef_state` : ReefState
- `ts` : Time step
- `loc` : Location
- `grp` : Functional group
- `h²` : Heritability
"""
function update_coral_tolerances!(
    reef_state::ReefState,
    ts::Int64,
    loc::Int64,
    grp::Int64,
    n_recruits::Int64;
    h²::Float32=0.3f0
)::Nothing
    ts2 = ts - 2
    if ts2 <= 0
        ts2 = 1
    end

    ts1 = ts - 1
    if ts1 <= 0
        ts1 = 1
    end

    # Get previous tolerance means
    wild_mean = reef_state.wild_dhw_tolerances[ts2, loc, grp, At(:mean)]
    deployed_mean = reef_state.deployed_dhw_tolerances[ts2, loc, grp, At(:mean)]

    # Calculate mixing proportions
    wild_prop, deployed_prop = calculate_inheritance_proportions(reef_state, ts2, loc, grp)

    # Calculate new mean based on mixing
    prev_mean_mixed = (wild_mean * wild_prop) + (deployed_mean * deployed_prop)

    # Get "current" tolerance means
    wild_mean = reef_state.wild_dhw_tolerances[ts1, loc, grp, At(:mean)]
    deployed_mean = reef_state.deployed_dhw_tolerances[ts1, loc, grp, At(:mean)]

    # Calculate mixing proportions
    wild_prop, deployed_prop = calculate_inheritance_proportions(reef_state, ts2, loc, grp)

    # Calculate new mean based on mixing
    mean_mixed = (wild_mean * wild_prop) + (deployed_mean * deployed_prop)

    # Apply breeder's equation to get new tolerance
    recruit_mean = breeders(prev_mean_mixed, mean_mixed, h²)

    # Weighted mean, based on current active population size
    pop = reef_state.wild_population[ts1, loc, grp]
    prop = n_recruits / (n_recruits + count(pop .> 0.0))
    new_grp_mean = Float32((recruit_mean * prop) + (prev_mean_mixed * (1.0 - prop)))

    # Update tolerances influenced by the new recruits
    return update_dhw_tol_mean!(reef_state, ts, loc, grp, new_grp_mean)
end

"""
    run_example!(reef_state::ReefState, env_conditions::YAXArray; rng::AbstractRNG=Random.GLOBAL_RNG)

Run example with constant incoming larvae.
"""
function run_example!(
    reef_state::ReefState,
    env_conditions::YAXArray;
    recruits=0.06f0,
    self_seed=0.3f0,
    rng::AbstractRNG=Random.GLOBAL_RNG
)
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
    # n_recruits = floor(Int, coral_sample_size * 0.2)
    # total_recruits = n_deploy + n_recruits
    recruit_μ = 1.5
    recruit_σ = 0.2

    # We want largest to smallest size classes for stratified sampling
    # (in `resample_population()`)
    class_diams = reverse.(CoralFlow.diameter_size_classes())

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

    # TODO: Refactor into `run_timestep()`
    for ts in timesteps[2:end]
        prev_ts = ts - 1

        # Used to determine available space, which affects settlement rate, etc.
        total_covers::Vector{Float32} = coral_cover(reef_state, prev_ts)

        recruits = fill(Float32[], n_locs, n_grps)  # Reset cache

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
                recruits[loc, grp] = Float32.(rand(rng, recruit_dist, n_loc_recruits))
                update_coral_tolerances!(reef_state, ts, loc, grp, n_loc_recruits)
            else
                # Have to update the previous time step's entry as later calculations
                # update tolerances influenced by the new recruits

                # Get "current" tolerance means
                wild_mean = reef_state.wild_dhw_tolerances[prev_ts, loc, grp, At(:mean)]
                deployed_mean = reef_state.deployed_dhw_tolerances[
                    prev_ts, loc, grp, At(:mean)
                ]

                reef_state.wild_dhw_tolerances[ts, loc, grp, At(:mean)] = wild_mean
                reef_state.deployed_dhw_tolerances[ts, loc, grp, At(:mean)] = deployed_mean
            end

            # TODO: Distribution should only be affected once juveniles reach maturity
            if reef_state.deployment_times[ts, loc, grp] > 0
                n_deploy = Int64(reef_state.deployment_times[ts, loc, grp])
                deploy_corals!(reef_state, ts, loc, n_deploy, grp)

                # Currently assuming the mean is not affected...
                # reef_state.deployed_dhw_tolerances[ts, loc, grp, At(:mean)] = val
            end
        end

        dhws = env_conditions[ts, :, At(:dhw)].data
        cyclone_cats = env_conditions[ts, :, At(:cyclone_category)].data

        for loc in 1:n_locs, grp in 1:n_grps
            pop_buffer .= 0.0f0  # Reset buffer

            fill_population_buffer!(
                reef_state, prev_ts, loc, grp, recruits[loc, grp], pop_buffer
            )

            # Diameters for entire population including recruits
            with_recruits = pop_buffer[pop_buffer .> 0.0f0]

            # Background mortality
            # TODO: apply location specific survival scaler
            apply_survival!(reef_state, grp, with_recruits)

            # Bleaching mortality
            new_mean, new_std, area_lost = bleaching_mortality!(
                with_recruits,
                dhws[loc],
                depth_coeffs[loc],
                reef_state.wild_dhw_tolerances[ts, loc, grp, :]
            )
            reef_state.wild_dhw_tolerances[ts, loc, grp, :] .= (new_mean, new_std)

            # Cyclone mortality
            # cyclone_probs = cyclone_mortality_prob.(p_sample, [cyclone_cat])
            # total_probs = cyclone_probs  # TODO: Add other probabilities
            # survivors_mask = rand(length(p_sample)) .> total_probs

            # Apply stratified sampling to maintain tracking of corals at expected
            # max density.
            # Now unnecessary as we're constraining to max density anyway
            # resample_wild_population!(
            #     reef_state, ts, loc, grp, with_recruits, class_diams[grp]
            # )
            update_pop_cache!(reef_state, with_recruits, loc)

            next_pop = @view(reef_state._pop_cache[loc, :])
            update_wild_sample!(reef_state, ts, loc, grp, next_pop[next_pop .> 0.0])
        end

        # Survivers grow...
        for grp in 1:n_grps
            # Wild population
            apply_growth!(
                reef_state,
                grp,
                inflection_points[grp],
                reef_state.wild_population[ts, :, grp],
                total_covers ./ reef_state.carrying_capacity
            )

            # Update total cover
            total_covers = coral_cover(reef_state, ts)

            # Deployed population
            apply_growth!(
                reef_state,
                grp,
                inflection_points[grp],
                reef_state.deployed_population[ts, :, grp],
                total_covers ./ reef_state.carrying_capacity
            )

            # Update total cover
            total_covers = coral_cover(reef_state, ts)
        end
    end

    return nothing
end

"""
    run_example(; n_ts=75, n_locs=100, with_dhw=true, growth_models=growth_models, survival_models=survival_models)

Run example with constant incoming larvae.
"""
function run_example(;
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
    example_env = generate_example_environment(n_ts, n_locs; with_dhw=with_dhw)

    run_example!(reef_state, example_env)

    return reef_state, example_env
end
