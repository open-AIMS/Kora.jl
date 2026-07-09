function apply_growth!(
    reef_state::ReefState, grp::Int64, inflection_point::F, diams::Vector{Vector{F}},
    reef_cover::Vector{F}
)::Nothing where {F<:Float32}
    model = reef_state.growth_models[grp]
    for i in eachindex(reef_cover)
        if isempty(diams[i])
            continue
        end
        constraint = space_constraint(reef_cover[i], 20.0f0, inflection_point)
        scaler = constraint * reef_state.location_scalers[1, i, grp]
        @inbounds for j in eachindex(diams[i])
            old_diam = diams[i][j]
            diams[i][j] = old_diam + (model(old_diam) * scaler)
        end
    end
    return nothing
end

"""
   apply_bleaching!(
       reef_state::ReefState,
       ts::Int64,
       loc::Int64,
       grp::Int64,
       dhw::Float32,
       depth_coeff::Float32
   )::Tuple{Float32,Float32}

Apply bleaching mortality to both wild and deployed coral populations.

# Returns
Tuple of (wild_area_lost, deployed_area_lost)
"""
function apply_bleaching!(
    reef_state::ReefState,
    ts::Int64,
    loc::Int64,
    grp::Int64,
    dhw::Float32,
    depth_coeff::Float32
)::Tuple{Float32,Float32}
    # Handle wild population
    wild_diams = wild_population(reef_state, ts - 1, loc, grp)
    wild_tols = @view(reef_state.wild_dhw_tolerances[ts - 1, loc, grp, :])
    wild_μ, wild_σ, wild_area_lost = bleaching_mortality!(
        wild_diams, dhw, depth_coeff, wild_tols, grp
    )

    # Update wild tolerances if there was mortality
    if wild_area_lost > 0.0f0
        reef_state.wild_dhw_tolerances[ts, loc, grp, 1] = wild_μ
        reef_state.wild_dhw_tolerances[ts, loc, grp, 2] = wild_σ
    end

    # Handle deployed population
    deployed_diams = deployed_population(reef_state, ts - 1, loc, grp)
    deployed_tols = @view(reef_state.deployed_dhw_tolerances[ts - 1, loc, grp, :])
    deployed_μ, deployed_σ, deployed_area_lost = bleaching_mortality!(
        deployed_diams, dhw, depth_coeff, deployed_tols, grp
    )

    # Update deployed tolerances if there was mortality
    if deployed_area_lost > 0.0f0
        reef_state.deployed_dhw_tolerances[ts, loc, grp, 1] = deployed_μ
        reef_state.deployed_dhw_tolerances[ts, loc, grp, 2] = deployed_σ
    end

    return wild_area_lost, deployed_area_lost
end

function apply_survival!(
    reef_state::ReefState, grp::Int64, diams::AbstractVector{F}, rng::AbstractRNG
)::Nothing where {F<:Float32}
    model = reef_state.survival_models[grp]
    survival!(model, diams, rng)

    return nothing
end

"""
    larval_production(
        reef_state::ReefState,
        ts::Int64,
        loc::Int64,
        grp::Int64
    )::Float32

Determine production of larvae for a population at a given location and time.
Only considers colonies above maturity threshold.
"""
function larval_production(
    reef_state::ReefState,
    maturity_thresholds::Vector{Float32},
    ts::Int64,
    loc::Int64,
    grp::Int64
)::Int64
    pop = coral_population(reef_state, ts, loc, grp)
    threshold = maturity_thresholds[grp]

    total = 0.0f0
    @simd for d in pop
        if d >= threshold
            total += larval_production(d, grp)
        end
    end
    return floor(Int64, total)
end

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
Only mature corals contribute to offspring tolerance inheritance.

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
    # Fast path: if no deployed corals exist, skip mixing calculations
    if all(iszero, reef_state.deployed_population)
        return _update_coral_tolerances_wild_only!(
            reef_state, ts, loc, grp, n_recruits, h²
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

    # Get t-2 tolerance means
    wild_mean_t2 = reef_state.wild_dhw_tolerances[ts2, loc, grp, 1]
    deployed_mean_t2 = reef_state.deployed_dhw_tolerances[ts2, loc, grp, 1]

    # Calculate mixing proportions at t-2
    wild_prop_t2, deployed_prop_t2 = calculate_inheritance_proportions(
        reef_state, ts2, loc, grp
    )
    prev_mean_mixed = (wild_mean_t2 * wild_prop_t2) + (deployed_mean_t2 * deployed_prop_t2)

    # Get t-1 tolerance means
    wild_mean_t1 = reef_state.wild_dhw_tolerances[ts1, loc, grp, 1]
    deployed_mean_t1 = reef_state.deployed_dhw_tolerances[ts1, loc, grp, 1]

    # Calculate mixing proportions at t-1
    wild_prop_t1, deployed_prop_t1 = calculate_inheritance_proportions(
        reef_state, ts1, loc, grp
    )
    mean_mixed = (wild_mean_t1 * wild_prop_t1) + (deployed_mean_t1 * deployed_prop_t1)

    # Apply breeder's equation
    recruit_mean = breeders(prev_mean_mixed, mean_mixed, h²)

    # Only mature corals reproduce; filter by maturity size
    mature_size = susceptibility_size_thresholds()[grp]
    wild_pop_t1 = reef_state.wild_population[ts1, loc, grp]
    deployed_pop_t1 = reef_state.deployed_population[ts1, loc, grp]
    n_existing_mature =
        count(wild_pop_t1 .>= mature_size) + count(deployed_pop_t1 .>= mature_size)

    # Weighted mean based on recruitment proportion
    prop = n_recruits / (n_recruits + n_existing_mature)
    new_grp_mean = Float32((recruit_mean * prop) + (mean_mixed * (1.0 - prop)))

    return update_dhw_tol_mean!(reef_state, ts, loc, grp, new_grp_mean)
end

function _update_coral_tolerances_wild_only!(
    reef_state::ReefState,
    ts::Int64,
    loc::Int64,
    grp::Int64,
    n_recruits::Int64,
    h²::Float32
)::Nothing
    ts2 = ts - 2
    if ts2 <= 0
        ts2 = 1
    end

    ts1 = ts - 1
    if ts1 <= 0
        ts1 = 1
    end

    # Get wild tolerance means only
    wild_mean_t2 = reef_state.wild_dhw_tolerances[ts2, loc, grp, 1]
    wild_mean_t1 = reef_state.wild_dhw_tolerances[ts1, loc, grp, 1]

    # Apply breeder's equation with wild population only
    recruit_mean = breeders(wild_mean_t2, wild_mean_t1, h²)

    # Only mature corals reproduce; filter by maturity size
    mature_size = susceptibility_size_thresholds()[grp]
    wild_pop_t1 = reef_state.wild_population[ts1, loc, grp]
    n_existing_mature = count(wild_pop_t1 .>= mature_size)

    # Weighted mean based on recruitment proportion
    prop = n_recruits / (n_recruits + n_existing_mature)
    new_grp_mean = Float32((recruit_mean * prop) + (wild_mean_t1 * (1.0 - prop)))

    return update_dhw_tol_mean!(reef_state, ts, loc, grp, new_grp_mean)
end
