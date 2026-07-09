# Positional-only overloads and WASM-specific entry helpers for WasmTarget compilation.
#
# WasmTarget emits wrong GC type indices when compiling Julia kwarg dispatch machinery
# (sym_in / _throw_argerror). These methods either add genuinely new positional signatures
# (different arg count from existing kwarg methods) or use renamed _wasm_* prefixes to
# avoid the method-overwrite restriction during Module precompilation.
#
# Included after all other reef source files so the originals are already defined.

# ── space_constraint (3-arg) ──────────────────────────────────────────────────
# Original: space_constraint(x, k; x0=0.96f0)  [2 positional + kwarg]
# New:      space_constraint(x, k, x0)          [3 positional — genuinely new]

@inline function space_constraint(x::F, k::F, x0::F)::F where {F<:Float32}
    if x >= 1.0f0
        return 0.0f0
    end
    return 1.0f0 / (1.0f0 + _euler_f32b^(k * (x - x0)))
end

# ── deploy_corals! (6-arg) ────────────────────────────────────────────────────
# Original: deploy_corals!(reef, ts, loc, n, grp; rng=GLOBAL_RNG)  [5 positional]
# New:      deploy_corals!(reef, ts, loc, n, grp, rng)              [6 positional]

function deploy_corals!(
    reef_state::ReefState, ts::Int64, loc::Int64, n::Int64, grp::Int64, rng::AbstractRNG
)::Nothing
    size_dist = size_distribution()[grp]
    edges = bin_edges()[grp, :]
    μ_g, σ_g = size_dist
    deploy_sample = _sample_lognormal_bounded(
        μ_g, σ_g, 0.0, Float64(maximum(edges)), n, rng
    )
    return update_deployed_sample!(reef_state, ts, loc, grp, deploy_sample)
end

# ── update_coral_tolerances! (6-arg) ─────────────────────────────────────────
# Original: update_coral_tolerances!(reef, ts, loc, grp, n; h²=0.3f0)  [5 positional]
# New:      update_coral_tolerances!(reef, ts, loc, grp, n, h_sq)       [6 positional]

function update_coral_tolerances!(
    reef_state::ReefState,
    ts::Int64,
    loc::Int64,
    grp::Int64,
    n_recruits::Int64,
    h_sq::Float32
)::Nothing
    has_deployed = false
    @inbounds for i in eachindex(reef_state.deployed_population)
        if !isempty(reef_state.deployed_population[i])
            has_deployed = true
            break
        end
    end
    if !has_deployed
        return _update_coral_tolerances_wild_only!(
            reef_state, ts, loc, grp, n_recruits, h_sq
        )
    end
    ts2 = ts - 2
    if ts2 <= 0
        ts2 = 1
    end
    ts1 = ts - 1
    if ts1 <= 0
        ts1 = 1
    end
    wild_mean_t2 = @inbounds reef_state.wild_dhw_tolerances[ts2, loc, grp, 1]
    deployed_mean_t2 = @inbounds reef_state.deployed_dhw_tolerances[ts2, loc, grp, 1]
    wild_prop_t2, deployed_prop_t2 = calculate_inheritance_proportions(
        reef_state, ts2, loc, grp
    )
    prev_mean_mixed = (wild_mean_t2 * wild_prop_t2) + (deployed_mean_t2 * deployed_prop_t2)
    wild_mean_t1 = @inbounds reef_state.wild_dhw_tolerances[ts1, loc, grp, 1]
    deployed_mean_t1 = @inbounds reef_state.deployed_dhw_tolerances[ts1, loc, grp, 1]
    wild_prop_t1, deployed_prop_t1 = calculate_inheritance_proportions(
        reef_state, ts1, loc, grp
    )
    mean_mixed = (wild_mean_t1 * wild_prop_t1) + (deployed_mean_t1 * deployed_prop_t1)
    recruit_mean = breeders(prev_mean_mixed, mean_mixed, h_sq)
    mature_size = susceptibility_size_thresholds()[grp]
    wild_pop_t1 = @inbounds reef_state.wild_population[ts1, loc, grp]
    deployed_pop_t1 = @inbounds reef_state.deployed_population[ts1, loc, grp]
    n_existing_mature =
        count(wild_pop_t1 .>= mature_size) + count(deployed_pop_t1 .>= mature_size)
    prop = n_recruits / (n_recruits + n_existing_mature)
    new_grp_mean = Float32((recruit_mean * prop) + (mean_mixed * (1.0 - prop)))
    return update_dhw_tol_mean!(reef_state, ts, loc, grp, new_grp_mean)
end

# ── initialize_coral_population! (4-arg) ─────────────────────────────────────
# Original: initialize_coral_population!(reef, loc, target; group_proportions, rng)  [3 positional]
# New:      initialize_coral_population!(reef, loc, target, rng)                     [4 positional]

function initialize_coral_population!(
    reef_state::ReefState, loc::Int64, target_pop_size::Int64, rng::AbstractRNG
)::Nothing
    group_props = Float32[0.1f0, 0.2f0, 0.25f0, 0.2f0, 0.25f0]
    edges = bin_edges()
    size_dist_all = size_distribution()
    for grp in 1:n_groups(reef_state)
        n_samples = round(Int, target_pop_size * group_props[grp])
        mu_g, sigma_g = size_dist_all[grp]
        initial_population = _sample_lognormal_bounded(
            mu_g, sigma_g, 0.0, Float64(maximum(edges[grp, :])), n_samples, rng
        )
        update_wild_sample!(reef_state, 1, loc, grp, initial_population)
    end
    n_locs2 = size(reef_state.wild_dhw_tolerances, 2)
    n_ts2 = size(reef_state.wild_dhw_tolerances, 1)
    for loc2 in 1:n_locs2
        reef_state.wild_dhw_tolerances[1, loc2, 1, 1] = 3.751612251
        reef_state.wild_dhw_tolerances[1, loc2, 2, 1] = 4.081622683
        reef_state.wild_dhw_tolerances[1, loc2, 3, 1] = 4.487465256
        reef_state.wild_dhw_tolerances[1, loc2, 4, 1] = 6.165751937
        reef_state.wild_dhw_tolerances[1, loc2, 5, 1] = 7.153507902
    end
    for ts in 1:n_ts2, loc2 in 1:n_locs2
        reef_state.wild_dhw_tolerances[ts, loc2, 1, 2] = 2.904433676
        reef_state.wild_dhw_tolerances[ts, loc2, 2, 2] = 3.159922076
        reef_state.wild_dhw_tolerances[ts, loc2, 3, 2] = 3.474118416
        reef_state.wild_dhw_tolerances[ts, loc2, 4, 2] = 4.773419097
        reef_state.wild_dhw_tolerances[ts, loc2, 5, 2] = 5.538122776
    end
    return nothing
end

# ── initialize_reef (6-arg positional) ───────────────────────────────────────
# Original: initialize_reef(; n_timesteps, n_locs, ...) [all-kwarg, 0 positional]
# New:      initialize_reef(n_ts, n_locs, area, density, gm, sm) [6 positional]

function initialize_reef(
    n_timesteps::Int64,
    n_locs::Int64,
    area::Float64,
    density::Int64,
    growth_models::AbstractCoralBehavior,
    survival_models::AbstractCoralBehavior
)::ReefState
    n_grps = length(TARGET_GROUPS)
    depths_val = 9.0f0

    location_scalers = fill(1.0f0, 2, n_locs, n_grps)
    wild_population = Array{Vector{Float32},3}(undef, n_timesteps, n_locs, n_grps)
    deployed_population = Array{Vector{Float32},3}(undef, n_timesteps, n_locs, n_grps)
    wild_population .= [Float32[]]
    deployed_population .= [Float32[]]
    deployment_times = zeros(Float32, n_timesteps, n_locs, n_grps)
    wild_dhw_tol = zeros(Float32, n_timesteps, n_locs, n_grps, 2)
    deployed_dhw_tol = zeros(Float32, n_timesteps, n_locs, n_grps, 2)
    mort = zeros(Float32, n_timesteps, n_locs, n_grps, 2)

    depths_vec = fill(depths_val, n_locs)
    area_vec = fill(Float32(area), n_locs)
    density_vec = fill(density, n_locs)

    buffer_size = Int64(maximum(area_vec) * maximum(density_vec) * 5)
    reef_state = ReefState(
        wild_population,
        deployed_population,
        deployment_times,
        growth_models.models,
        survival_models.models,
        location_scalers,
        density_vec,
        depths_vec,
        wild_dhw_tol,
        deployed_dhw_tol,
        mort,
        area_vec,
        10_000,
        zeros(Float32, n_locs, buffer_size),
        zeros(Float32, n_locs, buffer_size),
        zeros(Float32, n_locs),
        zeros(Float32, buffer_size),
        zeros(Float32, n_locs),
        zeros(Float32, n_locs)
    )
    return reef_state
end

# ── initialize_reef (7-arg positional, with depth) ───────────────────────────
# Original 6-arg: initialize_reef(n_ts, n_locs, area, density, gm, sm)
# New 7-arg:      initialize_reef(n_ts, n_locs, area, density, depth, gm, sm)

function initialize_reef(
    n_timesteps::Int64,
    n_locs::Int64,
    area::Float64,
    density::Int64,
    depth::Float64,
    growth_models::Kora.PolyGrowthModel{Float32},
    survival_models::Kora.PolySurvivalModel{Float32}
)::ReefState
    Kora.initialize_reef(;
        n_timesteps=Int(n_timesteps), n_locs=Int(n_locs),
        area=area, density=Int(density), depths=depth,
        growth_models=growth_models, survival_models=survival_models
    )
end

# ── run_model! (5-arg positional) ────────────────────────────────────────────
# Original: run_model!(reef, dhw; recruits, self_seed, rng)  [2 positional]
# New:      run_model!(reef, dhw, recruits, self_seed, rng)  [5 positional]
#
# Differences from the kwarg version:
#   - deploy_corals! called with 6 positional args (rng explicit)
#   - update_coral_tolerances! called with 6 positional args (h_sq=0.3f0 explicit)
#   - apply_growth! dispatches to the modified no-kwarg version in reef_dynamics.jl

function run_model!(
    reef_state::ReefState,
    dhw::Matrix{Float32},
    recruits_rate::Float32,
    self_seed_rate::Float32,
    rng::AbstractRNG
)::Nothing
    reset!(reef_state)

    timesteps::UnitRange{Int64} = 1:n_timesteps(reef_state)
    n_locs::Int64 = n_locations(reef_state)
    n_grps::Int64 = n_groups(reef_state)

    recruitment_proportion = recruits_rate
    self_seeding_proportion = self_seed_rate

    n_deploy = 0
    recruit_mu = 1.5f0
    recruit_sigma = 0.2f0

    carrying_cap = reef_state.carrying_capacity
    total_possible_colonies = carrying_cap .* reef_state.density

    inflection_points = growth_inflection_point()
    maturity_thresholds = mature_size_thresholds()

    recruits = fill(Float32[], n_locs, n_grps)
    depth_coeffs = depth_coefficient.(reef_state.depths)
    pop_buffer = reef_state._pop_buffer

    reset!(reef_state)

    @inbounds for loc in 1:n_locs, grp in 1:n_grps
        pop_buffer .= 0.0f0

        fill_population_buffer!(
            reef_state, 1, loc, grp, recruits[loc, grp], pop_buffer
        )

        with_recruits = pop_buffer[pop_buffer .> 0.0f0]

        apply_survival!(reef_state, grp, with_recruits, rng)

        tols = @inbounds Vector{Float32}([
            reef_state.wild_dhw_tolerances[1, loc, grp, 1],
            reef_state.wild_dhw_tolerances[1, loc, grp, 2]
        ])
        new_mean, new_std, area_lost = bleaching_mortality!(
            with_recruits,
            dhw[1, loc],
            depth_coeffs[loc],
            tols,
            grp
        )
        reef_state.wild_dhw_tolerances[1, loc, grp, 1] = new_mean
        reef_state.wild_dhw_tolerances[1, loc, grp, 2] = new_std

        update_pop_cache!(reef_state, with_recruits, loc)

        next_pop = reef_state._pop_cache[loc, :]
        update_wild_sample!(reef_state, 1, loc, grp, next_pop[next_pop .> 0.0])
    end

    recruits = fill(Float32[], n_locs, n_grps)

    for ts in timesteps[2:end]
        prev_ts = ts - 1

        total_covers::Vector{Float32} = coral_cover(reef_state, prev_ts)

        for i in eachindex(recruits)
            empty!(recruits[i])
        end

        @inbounds for loc in 1:n_locs, grp in 1:n_grps
            prod = larval_production(reef_state, maturity_thresholds, prev_ts, loc, grp)
            prod = prod * self_seeding_proportion * recruitment_proportion

            all_pop = total_population(reef_state, prev_ts, loc)
            avail_d_for_recruitment = all_pop < total_possible_colonies[loc]
            if avail_d_for_recruitment
                available_space = max(carrying_cap[loc] - total_covers[loc], 0.0)
                prod = (prod / reef_state.carrying_capacity[loc]) * available_space
                settlement_prop = min.(
                    (available_space * 50),
                    prod
                )
                n_loc_recruits = floor(Int64, settlement_prop)
            else
                n_loc_recruits = 0
            end

            if n_loc_recruits > 0
                recruits[loc, grp] = rand_truncated_normal(
                    rng, recruit_mu, recruit_sigma, 0.5f0, 2.5f0, n_loc_recruits
                )
                update_coral_tolerances!(reef_state, ts, loc, grp, n_loc_recruits, 0.3f0)
            else
                reef_state.wild_dhw_tolerances[ts, loc, grp, 1] = reef_state.wild_dhw_tolerances[
                    prev_ts, loc, grp, 1
                ]
                reef_state.deployed_dhw_tolerances[ts, loc, grp, 1] = reef_state.deployed_dhw_tolerances[
                    prev_ts, loc, grp, 1
                ]
            end

            if reef_state.deployment_times[ts, loc, grp] > 0
                n_deploy = Int64(reef_state.deployment_times[ts, loc, grp])
                deploy_corals!(reef_state, ts, loc, n_deploy, grp, rng)
            end
        end

        @inbounds for loc in 1:n_locs, grp in 1:n_grps
            pop_buffer .= 0.0f0

            fill_population_buffer!(
                reef_state, prev_ts, loc, grp, recruits[loc, grp], pop_buffer
            )

            with_recruits = pop_buffer[pop_buffer .> 0.0f0]

            apply_survival!(reef_state, grp, with_recruits, rng)

            tols = @inbounds Vector{Float32}([
                reef_state.wild_dhw_tolerances[ts, loc, grp, 1],
                reef_state.wild_dhw_tolerances[ts, loc, grp, 2]
            ])
            new_mean, new_std, area_lost = bleaching_mortality!(
                with_recruits,
                dhw[ts, loc],
                depth_coeffs[loc],
                tols,
                grp
            )
            reef_state.wild_dhw_tolerances[ts, loc, grp, 1] = new_mean
            reef_state.wild_dhw_tolerances[ts, loc, grp, 2] = new_std

            update_pop_cache!(reef_state, with_recruits, loc)

            next_pop = reef_state._pop_cache[loc, :]
            update_wild_sample!(reef_state, ts, loc, grp, next_pop[next_pop .> 0.0])
        end

        total_covers = coral_cover(reef_state, ts - 1)
        cover_proportions = total_covers ./ reef_state.carrying_capacity

        any_above = false
        @inbounds for v in cover_proportions
            if v > 0.85f0
                any_above = true
                break
            end
        end
        if any_above
            @inbounds for grp in 1:n_grps
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
            @inbounds for grp in 1:n_grps
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

        total_covers = coral_cover(reef_state, ts)
        @inbounds for loc in 1:n_locs
            if total_covers[loc] > reef_state.carrying_capacity[loc]
                scale = sqrt(reef_state.carrying_capacity[loc] / total_covers[loc])
                scale *= 1.0f0 - 8.0f0 * eps(Float32)
                @inbounds for grp in 1:n_grps
                    reef_state.wild_population[ts, loc, grp] .*= scale
                    reef_state.deployed_population[ts, loc, grp] .*= scale
                end
            end
        end
    end

    return nothing
end

# ── _wasm_generate_dhw ────────────────────────────────────────────────────────
# Renamed from generate_example_dhw(n, m) to avoid overwriting the kwarg version.
# Uses all-default parameters; replaces range(; length=) kwarg with inline math.

function _wasm_generate_dhw(n_years::Int64, n_locations::Int64)::Matrix{Float32}
    warming_rate::Float32 = 0.15f0
    seasonal_amplitude::Float32 = 1.2f0
    dhw_threshold::Float32 = 4.0f0
    noise_amplitude::Float32 = 0.9f0
    rng = Random.default_rng()
    n_years_f = Float32(n_years - 1)

    dhw_data = zeros(Float32, n_years, n_locations)

    @inbounds for loc in 1:n_locations
        spatial_offset = rand(rng) * 0.8f0
        weeks_above_threshold = 0

        @inbounds for t in 1:n_years
            yr_f = Float32(t - 1)
            warming_trend = warming_rate * yr_f
            seasonal_scaling = 1.0f0 + 0.2f0 * (yr_f / n_years_f)
            seasonal_cycle =
                seasonal_amplitude * seasonal_scaling *
                sin(Float32(2.0 * pi) * Float32(t % 12) / 12.0f0)
            weather_scaling = 1.0f0 + 0.3f0 * (yr_f / n_years_f)
            weather_noise = randn(rng) * noise_amplitude * weather_scaling
            temp_anomaly = (
                warming_trend + seasonal_cycle + weather_noise + spatial_offset
            )
            prev_dhw = t > 1 ? @inbounds(dhw_data[t - 1, loc]) : 0.0f0

            if temp_anomaly > dhw_threshold
                weeks_above_threshold += 1
                dhw_accumulation = (temp_anomaly - dhw_threshold) / 4.0f0
                if weeks_above_threshold >= 2
                    acute_probability = min(0.2f0 * (1.0f0 + warming_trend), 0.4f0)
                    if rand(rng) < acute_probability
                        time_factor = yr_f / n_years_f
                        base_spike = 3.0f0 + 2.0f0 * time_factor
                        spike_magnitude = rand(rng) * base_spike + 2.0f0
                        dhw_accumulation += spike_magnitude
                    end
                end
                base_dhw = prev_dhw * 0.92f0
                raw_dhw = base_dhw + dhw_accumulation
                soft_cap = 20.0f0
                if raw_dhw > soft_cap
                    excess = raw_dhw - soft_cap
                    damping_factor = 1.0f0 / (1.0f0 + 0.3f0 * excess)
                    time_factor = yr_f / n_years_f
                    max_fluctuation = 4.0f0 + 2.0f0 * time_factor
                    fluctuation =
                        (rand(rng) - 0.5f0) * min(max_fluctuation, excess * 0.6f0)
                    dhw_data[t, loc] = soft_cap + (excess * damping_factor) + fluctuation
                else
                    dhw_data[t, loc] = raw_dhw
                end
            else
                weeks_above_threshold = 0
                dhw_data[t, loc] = max(0.0f0, prev_dhw * 0.7f0)
            end
        end

        n_extreme_events = floor(Int, n_years / 12)
        @inbounds for _ in 1:n_extreme_events
            event_time = rand(rng, 1:n_years)
            time_progress = Float32(event_time - 1) / n_years_f
            event_probability = time_progress * 1.8f0
            if rand(rng, Float32) < event_probability
                base_duration = 2:5
                extra_duration = rand(rng, 0:floor(Int, 3 * time_progress))
                event_duration = rand(rng, base_duration) + extra_duration
                base_magnitude = 8.0f0 + 17.0f0 * time_progress
                event_magnitude = rand(rng) * 5.0f0 + base_magnitude
                @inbounds for t in event_time:min(event_time + event_duration, n_years)
                    relative_pos = (t - event_time) / event_duration
                    if relative_pos <= 0.3f0
                        scaling = relative_pos / 0.3f0
                    elseif relative_pos >= 0.7f0
                        scaling = 1.0f0 - ((relative_pos - 0.7f0) / 0.3f0)
                    else
                        scaling = 1.0f0
                    end
                    dhw_value = event_magnitude * Float32(scaling)
                    dhw_data[t, loc] = max(dhw_data[t, loc], dhw_value)
                end
            end
        end
    end

    return dhw_data
end

# ── _wasm_init_coral_pop! ─────────────────────────────────────────────────────
# Renamed from initialize_coral_population!(reef) to avoid overwriting kwarg version.

function _wasm_init_coral_pop!(reef_state::ReefState)::Nothing
    n_locs = n_locations(reef_state)
    sample_size = ceil(Int64, maximum(reef_state.carrying_capacity) * 5)
    rng = Random.default_rng()
    for loc in 1:n_locs
        initialize_coral_population!(reef_state, loc, sample_size, rng)
    end
    return nothing
end

# ── _wasm_run_ensemble! ───────────────────────────────────────────────────────
# Renamed from run_ensemble!(reef, dhw, params) to avoid overwriting kwarg version.
# Uses 5-arg positional run_model! with baked-in defaults.

function _wasm_run_ensemble!(
    reef_state::ReefState,
    dhw::Matrix{Float32},
    ensemble_params::Matrix{Float64}
)
    n_ensemble = size(ensemble_params, 2)
    n_ts = n_timesteps(reef_state)
    n_locs = n_locations(reef_state)
    n_grps = n_groups(reef_state)
    rng = Random.default_rng()

    ec, egc, ejc, ewdt = _alloc_ensemble_results(n_ts, n_locs, n_grps, n_ensemble)

    for i in 1:n_ensemble
        params_i = ensemble_params[:, i]
        set_population!(reef_state, params_i)
        run_model!(reef_state, dhw, 0.06f0, 0.3f0, rng)
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
