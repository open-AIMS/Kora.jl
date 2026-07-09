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
        size_dist = (
            (_x[1], _x[6]),
            (_x[2], _x[7]),
            (_x[3], _x[8]),
            (_x[4], _x[9]),
            (_x[5], _x[10])
        )
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
    mean_colony_cover_m2()::Float32

Expected cover (m²) per colony, averaged across all functional groups. Computed
analytically from the truncated log-normal size distributions used by
`initialize_coral_population!`. Useful for converting a target cover fraction
to an initial colony count.

For X ~ TruncLN(μ, σ, 0, b): E[X²] = exp(2μ + 2σ²) · Φ((ln b − μ − 2σ²)/σ) / Φ((ln b − μ)/σ),
where Φ is the standard normal CDF approximated via `rational_erf`.
"""
function mean_colony_cover_m2()::Float32
    dists = size_distribution()
    edges = bin_edges()
    n_grps = size(edges, 1)
    total = 0.0
    for grp in 1:n_grps
        μ, σ = Float64.(dists[grp])
        b = Float64(maximum(edges[grp, :]))
        lnb = log(b)
        σ² = σ * σ
        # Φ(x) = (1 + erf(x / √2)) / 2, using the bespoke rational_erf
        c = 0.7071067811865476
        phi_num = (1.0 + rational_erf(((lnb - μ - 2.0 * σ²) / σ) * c)) * 0.5
        phi_den = (1.0 + rational_erf(((lnb - μ) / σ) * c)) * 0.5
        ex2 = exp(2.0 * μ + 2.0 * σ²) * phi_num / phi_den
        total += ex2 * Float64(_pi_f32) * 0.25 * Float64(_cover_scale)
    end
    return Float32(total / n_grps)
end
