"""Methods to support calibration."""

"""
    assign_scalers!(reef_state::ReefState, x::Vector)::Nothing

Assign growth and survival scalers, assuming they are identical for each location.
"""
function assign_scalers!(reef_state::ReefState, x::Vector)::Nothing
    n_grps = n_groups(reef_state)

    # Modify location scalers
    growth_scalers = x[1:n_grps]
    # survival_scalers = x[(n_grps + 1):end]
    # collated = collect(zip(growth_scalers, survival_scalers))

    for i in 1:n_grps
        # reef_state.location_scalers[:, :, i] .= collated[i]
        reef_state.location_scalers[1, :, i] .= growth_scalers[i]
    end

    return nothing
end

"""
    set_population!(reef_state::ReefState, x::Vector)::Nothing

Set the initial population state.
"""
function set_population!(reef_state::ReefState, x::Vector)::Nothing
    pop_density = x[1]
    total_initial_pop = ceil(
        Int64, floor(Int64, pop_density) * reef_state.carrying_capacity[1]
    )

    if length(x) > 6
        # Get population size distribution
        _x = Float32.(x[7:16])
        size_dist = [
            LogNormal{Float32}(_x[1], _x[6]),  # tabular Acropora
            LogNormal{Float32}(_x[2], _x[7]),  # corymbose Acropora
            LogNormal{Float32}(_x[3], _x[8]),  # Pocillopora + non-Acropora corymbose
            LogNormal{Float32}(_x[4], _x[9]),   # Small massives and encrusting
            LogNormal{Float32}(_x[5], _x[10])    # Large massives
        ]
    else
        size_dist = Kora.size_distribution()
    end

    # Proportion of each functional group
    proportions = Float32.(x[2:6])
    initialize_coral_population!(
        reef_state, 1, total_initial_pop; group_proportions=proportions, size_dist=size_dist
    )

    return nothing
end

"""
    mean_colony_cover_m2(n_per_grp::Int=20_000)::Float32

Expected cover (m²) per colony, averaged across all functional groups, sampled
from the same truncated log-normal size distributions used by
`initialize_coral_population!`. Useful for converting a target cover fraction
to an initial colony count.
"""
function mean_colony_cover_m2(n_per_grp::Int=20_000)::Float32
    rng    = Random.MersenneTwister(0)
    dists  = size_distribution()
    edges  = bin_edges()
    n_grps = size(edges, 1)
    total  = 0.0f0
    for grp in 1:n_grps
        d      = truncated(dists[grp], 0.0f0, maximum(edges[grp, :]))
        samples = Float32.(rand(rng, d, n_per_grp))
        total  += sum(cover_cm_to_m2.(samples))
    end
    return total / (n_per_grp * n_grps)
end
