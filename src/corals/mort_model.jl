using CurveFit

# Private consts - only used in this file
const _euler_f32 = Float32(ℯ)

"""
    bleaching_susceptibility(x, k; x0=150.0)

Modify susceptibility of corals to bleaching as a function of its diameter.

# Arguments
- `x` : Proportion of available space currently taken
- `k` : Steepness parameter, where the larger the value, the faster the drop (suggest k=0.15).
- `x0` : Inflection point (based on diameter) where bleaching begins to have much less mortality (default: 150cm diameter).
"""
function bleaching_susceptibility(x::F, k::F; x0::F=150.0f0)::F where {F<:Float32}
    return 1.0f0 / (1.0f0 + _euler_f32^(k*(x - x0)))
end
function bleaching_susceptibility!(x::AbstractArray{F}, k::F, cache::AbstractArray{F}; x0::F=150.0f0)::Nothing where {F<:Float32}
    Threads.@threads for i in eachindex(x, cache)
        @inbounds cache[i] = bleaching_susceptibility(x[i], k; x0=x0)
    end

    return nothing
end

"""
    bleaching_mortality!(
        pop::AbstractMatrix{F},  # Population size distribution
        dhw::Vector{F},
        depth_coeff::Vector{F},
        tols::YAXArray{F,2},
        susceptibility_cache::AbstractMatrix{F},
        pop_cache::AbstractMatrix{F}
    ) where {F<:Float32}

Bleaching mortality.

# Arguments
- `pop` : Current sample population
- `dhw` : DHW for all locations
- `depth_coeff` : Depth coefficient for all locations
- `tols` : Current mean DHW tolerance
- `susceptibility_cache` : Cache to store intermediate susceptiblity values
- `pop_cache` : Cache to store intermediate population values

# Returns
Tuple: Updated mean DHW tolerance, cover lost (in m²)
"""
function bleaching_mortality!(
    pop::AbstractMatrix{F},  # Population size distribution
    dhw::Vector{F},
    depth_coeff::Vector{F},
    tols::YAXArray{F,2},
    susceptibility_cache::AbstractMatrix{F},
    pop_cache::AbstractMatrix{F}
)::Tuple where {F<:Float32}
    # Calculate adaptation-based population effect (from original bleaching_mortality)
    μ = tols[:, At(:mean)].data
    stdev = tols[:, At(:stdev)].data
    affected_prop::Vector{F} = truncated_normal_cdf.(
        dhw, μ, stdev,
        4.0f0,
        μ .+ 10.0f0
    )

    # Apply depth coefficient to affected proportion
    base_affected = affected_prop .* depth_coeff

    # Apply size-dependent mortality to each size class
    pop_cache .= pop

    # Calculate size-specific mortality modifier
    bleaching_susceptibility!(pop, 0.15f0, susceptibility_cache)

    # Explicit loop to avoid temporary allocations
    Threads.@threads for i in axes(pop, 2)
        @inbounds @simd for j in axes(pop, 1)
            if pop[j, i] > 5.0f0
                pop[j, i] *= (1.0f0 - (base_affected[j] * susceptibility_cache[j, i]))
            end
        end
    end

    cover_cm_to_m2!(pop_cache .- pop, pop_cache)
    area_lost = sum(pop_cache; dims=2)

    # Calculate new mean thermal tolerance if mortality occurred
    if any(area_lost .> 0.0f0)
        mort_locs = area_lost .> 0.0f0
        # Re-create distribution maintaining genetic variance
        μ[mort_locs] .= truncated_normal_mean.(
            μ[mort_locs], stdev[mort_locs], 4.0f0, μ[mort_locs] .+ 10.0f0
        )
    end

    return μ, area_lost
end

"""
    survival(model::Polynomial{F}, ext::F)::F where {F<:Float32}

Proportion of survival from background mortality.

# Arguments
- `model` : Polynomial
- `ext` : Current linear extension

# Returns
Survival from background mortality
"""
@inline function survival(model::Polynomial{F}, ext::F)::F where {F<:Float32}
    return clamp(model(ext), 0.0f0, 1.0f0)
end
function survival!(model::Polynomial{F}, ext::AbstractMatrix{F}, cache::AbstractMatrix{F})::Nothing where {F<:Float32}
    Threads.@threads for i in eachindex(ext, cache)
        @inbounds cache[i] = survival(model, ext[i])
    end

    return nothing
end

survival_models = Polynomial[]
model_coefs = []
rmse_scores = []
r2_scores = []
for (xi, yi) in zip(eachrow(bin_edges()), eachrow(survival_rates()))
    m = curve_fit(Polynomial, xi, yi, 3)
    push!(rmse_scores, RMSE(m.(xi), yi))
    push!(r2_scores, R2(m.(xi), vec(yi)))
    push!(survival_models, m)
    push!(model_coefs, m.coeffs)
end
