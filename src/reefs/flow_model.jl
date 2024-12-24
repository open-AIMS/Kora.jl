"""
    run_example!(reef_state::ReefState, env_conditions::YAXArray; rng::AbstractRNG=Random.GLOBAL_RNG)

Run example with constant incoming larvae.
"""
function run_example!(
    reef_state::ReefState, env_conditions::YAXArray; rng::AbstractRNG=Random.GLOBAL_RNG
)
    timesteps::UnitRange{Int64} = 1:n_timesteps(reef_state)
    n_locs::Int64 = n_locations(reef_state)
    n_grps::Int64 = n_groups(reef_state)
    coral_sample_size::Int64 = pop_sample_size(reef_state)
    set_sample_size::UnitRange{Int64} = 1:coral_sample_size

    # Deployment per group
    # Assuming a 200m study area, this comes to ~7 deployments per m²
    # which is relatively high deployment density
    # (n * n_grps) / area
    # (280 * 5) / 200 = 7
    n_deploy = 280
    n_recruits = floor(Int, coral_sample_size * 0.2)
    total_recruits = n_deploy + n_recruits
    recruit_μ = 1.0
    recruit_σ = 0.2

    _total_pop_sample_size::Int64 = total_recruits + coral_sample_size

    # We want largest to smallest size classes for stratified sampling
    # (in `update_size_distribution()`)
    class_diams = reverse.(CoralFlow.diameter_size_classes())

    # Mock recruitment
    # TODO: Remove once ready
    recruit_dist = truncated(Normal(recruit_μ, recruit_σ), 0.5, 2.5)

    # Convenience var
    carrying_cap = reef_state.carrying_capacity

    # Group modifiers that adjusts when growth slows down as it approaches a
    # proportionate limit relative to available space
    grp_mod = growth_inflection_point()

    cover_t::Vector{Float32} = coral_cover(reef_state, 1)
    recruits = zeros(Float32, n_locs, n_grps, total_recruits)
    depth_coeffs = depth_coefficient.(reef_state.depths)
    for ts in timesteps[2:end]
        prev_ts = ts - 1

        # Determine settlement rate based on available space
        # total_cover = coral_cover(reef_state, prev_ts, loc)
        total_covers::Vector{Float32} = cover_t

        # Some recruitment happens
        # Calculate across all locations and groups for time `ts`...
        for loc in 1:n_locs
            for grp in 1:n_grps
                recruits[loc, grp, :] .= 0.0  # Reset cache

                # TODO: Fix - massives tend to produce very large numbers of recruits
                # even when the current population is tiny. Could assume these are
                # external larvae, but that narrative is not internally consistent...
                # Scale recruitment by cover proportion
                prod = larval_production(reef_state, prev_ts, loc, grp)
                n_loc_recruits = min(
                    prod,
                    n_recruits
                )

                # Always deploy
                n_loc_recruits += n_deploy

                if n_loc_recruits > 0
                    recruits[loc, grp, 1:n_loc_recruits] .= rand(rng, recruit_dist, n_loc_recruits)
                end

                # TODO: Clean up - this handles DHW tolerance inheritance
                curr_dhw_tol = reef_state.dhw_tolerances[ts-1, loc, grp, 1]
                if (ts - 2) <= 0
                    prev_dhw_tol = reef_state.dhw_tolerances[ts-1, loc, grp, 1]
                else
                    if !(any(reef_state.dhw_tolerances[1:ts-1, loc, grp, 1].data .!= curr_dhw_tol))
                        prev_dhw_tol = reef_state.dhw_tolerances[ts-2, loc, grp, 1]
                    else
                        idx = last(findall(reef_state.dhw_tolerances[1:ts-2, loc, grp, 1].data .!= curr_dhw_tol))
                        prev_dhw_tol = reef_state.dhw_tolerances[idx, loc, grp, 1]
                    end
                end

                if curr_dhw_tol != prev_dhw_tol
                    rec_mean = breeders(prev_dhw_tol, curr_dhw_tol, 0.3f0)
                    pop = reef_state.pop_sample[prev_ts, 1, 1, :].data

                    # Weighted mean, based on current active population size
                    prop = n_loc_recruits / (n_loc_recruits + count(pop .> 0.0))
                    new_grp_mean = Float32((rec_mean * prop) + (prev_dhw_tol * (1.0 - prop)))
                else
                    new_grp_mean = curr_dhw_tol
                end

                # Have to update the previous time step's entry as later calculations
                # assume that's the input value.
                update_dhw_tol_mean!(reef_state, prev_ts, loc, grp, new_grp_mean)
            end
        end

        rec_cover::Vector{Float32} = recruit_cover(reef_state, recruits)
        exceed_capacity = (total_covers .+ rec_cover) .> carrying_cap
        if any(exceed_capacity)
            scale_factor = max.(carrying_cap[exceed_capacity] .- total_covers[exceed_capacity], 0.0f0) ./ carrying_cap[exceed_capacity]
            _n_recruits = max.(0, floor.(Int32, n_recruits * scale_factor))

            if any(_n_recruits .> 0)
                max_rec = [1:n for n in _n_recruits]
                for (n, max_n) in enumerate(max_rec)
                    recruits[exceed_capacity, :, max_n] .= recruits[exceed_capacity, :, max_n]
                    recruits[exceed_capacity, :, (_n_recruits[n]+1):end] .= 0.0f0
                end
            end
        end

        dhws = env_conditions[ts, :, At(:dhw)].data
        cyclone_cats = env_conditions[ts, :, At(:cyclone_category)].data

        curr_pop_cache = similar(reef_state._pop_buffer[:, 1:_total_pop_sample_size])
        pop_cache = similar(reef_state._pop_buffer[:, 1:_total_pop_sample_size])
        @inbounds for grp in 1:n_grps
            reef_state._pop_buffer .= 0.0f0  # Reset cache
            reef_state._pop_buffer[:, set_sample_size] .= reef_state.pop_sample.data[prev_ts, :, grp, :]
            reef_state._pop_buffer[:, (coral_sample_size+1):_total_pop_sample_size] .= recruits[:, grp, :]
            with_recruits = @views(reef_state._pop_buffer[:, 1:_total_pop_sample_size])

            # Apply mortality
            surv_model = reef_state.survival_models[grp]

            # Background mortality
            # (TODO: apply location specific survival scaler)
            survival!(surv_model, with_recruits, curr_pop_cache)

            # Bleaching mortality
            dhw_tols = @view(reef_state.dhw_tolerances[prev_ts, :, grp, :])
            curr_mean, area_lost = bleaching_mortality!(
                with_recruits, dhws, depth_coeffs, dhw_tols, curr_pop_cache, pop_cache
            )
            update_dhw_tol_mean!(reef_state, ts, grp, curr_mean)

            # Cyclone mortality
            # cyclone_probs = cyclone_mortality_prob.(p_sample, [cyclone_cat])
            # total_probs = cyclone_probs  # TODO: Add other probabilities
            # survivors_mask = rand(length(p_sample)) .> total_probs

            # update_mortalities!(reef_state, ts, loc, Float32[prop_mort, mean(total_probs)])

            scalers = reef_state.location_scalers[1, :, grp].data
            if grp == 1
                scalers .+= 0.0
            end
            growth!(
                reef_state.growth_models[grp],
                with_recruits,
                total_covers ./ reef_state.carrying_capacity,
                grp_mod[grp],
                scalers
            )

            # Apply stratified sampling to maintain desired sample size for all locations
            for loc in axes(with_recruits, 1)
                update_size_distribution!(reef_state, ts, loc, grp, with_recruits, class_diams[grp])
            end
        end

        cover_t = coral_cover(reef_state, ts)
    end

    return nothing
end

"""
    run_example(; n_ts=75, n_locs=100, with_dhw=true, growth_models=growth_models, survival_models=survival_models)

Run example with constant incoming larvae.
"""
function run_example(; n_ts=75, n_locs=100, with_dhw=true, area=100.0, pop_density=15.0, growth_models=growth_models, survival_models=survival_models)
    sample_size = ceil(Int64, area * pop_density) * 2
    reef_state = initialize_reef(; n_timesteps=n_ts, n_locs=n_locs, area=area, sample_size=sample_size, growth_models=growth_models, survival_models=survival_models)
    initialize_coral_population!(reef_state)
    example_env = generate_example_environment(n_ts, n_locs; with_dhw=with_dhw)

    run_example!(reef_state, example_env)

    return reef_state, example_env
end
