function apply_growth!(
    reef_state::ReefState, grp::Int64, inflection_point::F, diams::Vector{Vector{F}},
    reef_cover::Vector{F}
)::Nothing where {F<:Float32}
    scalers = reef_state.location_scalers[At(:growth), :, grp].data
    growth!(reef_state.growth_models[grp], diams, reef_cover, inflection_point, scalers)

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
        wild_diams, dhw, depth_coeff, wild_tols
    )

    # Update wild tolerances if there was mortality
    if wild_area_lost > 0.0f0
        reef_state.wild_dhw_tolerances[ts, loc, grp, At(:mean)] = wild_μ
        reef_state.wild_dhw_tolerances[ts, loc, grp, At(:stdev)] = wild_σ
    end

    # Handle deployed population
    deployed_diams = deployed_population(reef_state, ts - 1, loc, grp)
    deployed_tols = @view(reef_state.deployed_dhw_tolerances[ts - 1, loc, grp, :])
    deployed_μ, deployed_σ, deployed_area_lost = bleaching_mortality!(
        deployed_diams, dhw, depth_coeff, deployed_tols
    )

    # Update deployed tolerances if there was mortality
    if deployed_area_lost > 0.0f0
        reef_state.deployed_dhw_tolerances[ts, loc, grp, At(:mean)] = deployed_μ
        reef_state.deployed_dhw_tolerances[ts, loc, grp, At(:stdev)] = deployed_σ
    end

    return wild_area_lost, deployed_area_lost
end

function apply_survival!(
    reef_state::ReefState, grp::Int64, diams::AbstractVector{F}
)::Nothing where {F<:Float32}
    model = reef_state.survival_models[grp]
    survival!(model, diams)

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

    return floor(Int64, sum([larval_production(d, grp) for d in pop if d >= threshold]))
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
