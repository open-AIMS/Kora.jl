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
    depth_coefficient(d::Union{Int64,Float64})::Float64

Model by Baird et al., [1] providing an indication of a relationship between bleaching
and depth.

# Arguments
- `d` : median depth of location in meters

# Returns
Proportion of population affected by bleaching at depth `x`.
Values are constrained such that 0.0 <= x <= 1.0

# References
1. Baird, A., Madin, J., Álvarez-Noriega, M., Fontoura, L., Kerry, J., Kuo, C.,
     Precoda, K., Torres-Pulliza, D., Woods, R., Zawada, K., & Hughes, T. (2018).
   A decline in bleaching suggests that depth can provide a refuge from global
     warming in most coral taxa.
   Marine Ecology Progress Series, 603, 257-264.
   https://doi.org/10.3354/meps12732
"""
function depth_coefficient(d::F)::Float32 where {F<:Real}
    return Float32(max(min(exp(-0.07551 * (d - 2.0)), 1.0), 0.0))
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

Diameter size reduction due to bleaching mortality and partial mortality.

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
    if all(dhw .== 0.0)
        return tols[:, At(:mean)].data, 0.0f0
    end

    # Calculate adaptation-based population effect (from original bleaching_mortality)
    μ = tols[:, At(:mean)].data
    stdev = tols[:, At(:stdev)].data
    affected_prop::Vector{F} = truncated_normal_cdf.(
        dhw, μ, stdev,
        2.0f0,
        μ .+ 10.0f0
    )

    # Apply depth coefficient to affected proportion
    base_affected = affected_prop .* depth_coeff

    # Apply size-dependent mortality to each size class
    pop_cache .= pop

    # Calculate size-specific mortality modifier
    bleaching_susceptibility!(pop, 0.15f0, susceptibility_cache)

    # Reduction in size due to partial mortality or mortality
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
        # Re-create distribution, maintaining genetic variance
        μ[mort_locs] .= truncated_normal_mean.(
            μ[mort_locs], stdev[mort_locs], 2.0f0, μ[mort_locs] .+ 10.0f0
        )
    end

    return μ, area_lost
end

"""
    survival(model::M, diam::F)::F where {M,F<:Float32}

Proportion of survival from background mortality.

# Arguments
- `model` : Function
- `diam` : Current diameter

# Returns
Survival from background mortality
"""
@inline function survival(model::M, diam::F)::Bool where {M,F<:Float32}
    return rand(Float32) < model(diam)
end
function survival!(model::M, diam::AbstractMatrix{F}, cache::AbstractMatrix{F})::Nothing where {M,F<:Float32}
    Threads.@threads for i in eachindex(diam)
        @inbounds cache[i] = survival(model, diam[i])
    end

    diam .*= cache

    return nothing
end

survival_models = Function[]
mort_model_coefs = []
mort_rmse_scores = []
mort_r2_scores = []
for (xi, yi) in zip(eachrow(bin_edges()[:, 1:end-1]), eachrow(survival_rates()))
    xi, yi = log.(xi), log.(yi)
    m = curve_fit(Polynomial, xi, yi, 3)
    model = x -> x .< exp(xi[1]) ? exp(yi[1]) : clamp(exp(m(log(x))), 0.0f0, 1.0f0)
    push!(mort_rmse_scores, RMSE(model.(xi), yi))
    push!(mort_r2_scores, R2(model.(xi), vec(yi)))
    push!(survival_models, model)
    push!(mort_model_coefs, m.coeffs)
end

"""
    SurvivalModel <: AbstractCoral

Functional relationships between coral size and survival.
"""
@kwdef struct SurvivalModel <: AbstractCoral
    models::Vector{Function}
    rmse_scores::Vector{Float64} = []
    r2_scores::Vector{Float64} = []
end

function SurvivalModel(models)
    return SurvivalModel(models, Float32[], Float32[])
end
