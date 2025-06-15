function apply_growth!(
    reef_state::ReefState, grp::Int64, inflection_point::F, diams::Vector{Vector{F}}, reef_cover::Vector{F}
)::Nothing where {F<:Float32}
    scalers = reef_state.location_scalers[1, :, grp].data
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
    wild_tols = @view(reef_state.wild_dhw_tolerances[ts-1, loc, grp, :])
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
    deployed_tols = @view(reef_state.deployed_dhw_tolerances[ts-1, loc, grp, :])
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

function apply_survival!(reef_state::ReefState, grp::Int64, diams::AbstractVector{F})::Nothing where {F<:Float32}
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
function larval_production(reef_state::ReefState, maturity_thresholds::Vector{Float32}, ts::Int64, loc::Int64, grp::Int64)::Int64
    pop = coral_population(reef_state, ts, loc, grp)
    threshold = maturity_thresholds[grp]

    return floor(Int64, sum([larval_production(d, grp) for d in pop if d >= threshold]))
end