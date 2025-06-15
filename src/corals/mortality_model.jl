# Private consts - only used in this file
const _euler_f32 = Float32(ℯ)

"""
    bleaching_susceptibility(x, k; x0=150.0)

Modify susceptibility of corals to bleaching as a function of its diameter.

# Arguments
- `x` : Diameter of coral
- `k` : Steepness parameter, where the larger the value, the faster the drop (suggest k=0.15).
- `x0` : Inflection point (based on diameter) where bleaching begins to have much less mortality (default: 150cm diameter).
"""
function bleaching_susceptibility(x::F, k::F; x0::F=150.0f0)::F where {F<:Float32}
    return 1.0f0 / (1.0f0 + _euler_f32^(k * (x - x0)))
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
        diams::Vector{F},
        dhw::F,
        depth_coeff::F,
        tols::YAXArray{F, 1}
    )::Tuple where {F<:Float32}

Diameter size reduction due to bleaching mortality and partial mortality.

# Arguments
- `diams` : population diameter(s)
- `dhw` : Degree heating week experienced
- `depth_coeff` : Depth coefficient that ameliorates heat stress
- `tols` : population heat tolerances (mean and stdev)

# Returns
Tuple: Updated mean DHW tolerance, cover lost (in m²)
"""
function bleaching_mortality!(
    diams::AbstractVector{F},
    dhw::F,
    depth_coeff::F,
    tols::YAXArray{F,1}
)::Tuple where {F<:Float32}
    if all(dhw .< 4.0)
        return tols[At(:mean)], tols[At(:stdev)], 0.0f0
    end

    # Calculate adaptation-based population effect (from original bleaching_mortality)
    μ::F = tols[At(:mean)]
    stdev::F = tols[At(:stdev)]
    affected_prop::F = truncated_normal_cdf(
        dhw, μ, stdev,
        4.0f0,  # bleaching doesn't really occur until 4 DHW
        μ + 10.0f0
    )

    # Apply depth coefficient to affected proportion
    base_affected = affected_prop .* depth_coeff
    if base_affected == 0.0f0
        return μ, stdev, 0.0
    end

    diam_cache = copy(diams)

    # Apply size-dependent mortality to each size class
    # (the reduction in size due to partial mortality or mortality).
    # Using an explicit loop here to avoid temporary allocations
    Threads.@threads for i in eachindex(diams)
        if diams[i] > 2.5f0
            # Calculate size-specific mortality modifier
            diams[i] *= (1.0f0 - (base_affected * bleaching_susceptibility(diams[i], 0.15f0)))
        end
    end

    current_cover = cover_cm_to_m2(diams)
    cover_cm_to_m2!(max.(diam_cache .- diams, 0.0f0), diam_cache)
    area_lost = min(current_cover, sum(diam_cache))

    if any(area_lost > 0.0f0) && ((current_cover - area_lost) > 0.0f0)
        μ = truncated_normal_mean(μ, stdev, 4.0f0, μ + 10.0f0)
    end

    return μ, stdev, area_lost
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
function survival!(model::M, diam::AbstractVector{F})::Nothing where {M,F<:Float32}
    Threads.@threads for i in eachindex(diam)
        @inbounds diam[i] *= survival(model, diam[i])
    end

    return nothing
end

function log_likelihood(y_true, y_pred)
    return sum(y_true .* log.(y_pred .+ eps(Float32)) .+
               (1 .- y_true) .* log.(1 .- y_pred .+ eps(Float32)))
end

function null_log_likelihood(y_true)
    p_null = mean(y_true)
    return sum(y_true .* log(p_null) .+
               clamp.(1 .- y_true, 0.0, 1.0) .* log(clamp.(1 .- p_null, 0.0, 1.0)))
end

function mcfadden_r2(y_true, y_pred)
    return 1.0 - (log_likelihood(y_true, y_pred) / null_log_likelihood(y_true))
end

function brier_score(y_true, y_pred)
    return mean((y_pred .- y_true) .^ 2)
end

"""
    LogisticSurvivalModel <: AbstractCoralBehavior

Functional relationships between coral size and survival.
"""
struct LogisticSurvivalModel <: AbstractCoralBehavior
    models::Vector{Function}
    rmse_scores::Vector{Float32}
    log_likelihood_scores::Vector{Float32}
    mcfadden_r2_scores::Vector{Float32}
    brier_scores::Vector{Float32}

    function LogisticSurvivalModel()
        return new(
            survival_models,
            Float32.(mort_rmse_scores),
            Float32.(mort_ll_scores),
            Float32.(mort_mcfadden_scores),
            Float32.(mort_brier_scores)
        )
    end
end

function Base.show(io::IO, ::MIME"text/plain", x::LogisticSurvivalModel)
    println(io, "\nSurvival Model Performance Metrics:")
    println(io, "─"^40)

    group_names = ["Tabular Acropora", "Corymbose Acropora",
        "Corymbose non-Acropora", "Small massives", "Large massives"]

    explainer = """
        McFadden's R²: -∞ - 1.0; Higher is better.
        Log Likelihood: -∞ - 0.0; Higher is better.
        Brier Score: 0.0 - 1.0; Lower is better; 0.25 = random guess.
        """

    println(io, explainer)

    for i in eachindex(x.models)
        println(io, "\nGroup: $(group_names[i])")
        println(io, "RMSE:     $(Printf.@sprintf("%.3f", x.rmse_scores[i]))")
        println(io, "McFadden's R²:     $(Printf.@sprintf("%.3f", x.mcfadden_r2_scores[i]))")
        println(io, "Log Likelihood:     $(Printf.@sprintf("%.3f", x.log_likelihood_scores[i]))")
        println(io, "Brier Score:       $(Printf.@sprintf("%.3f", x.brier_scores[i]))")

        # println(io, "\nConfusion Matrix (threshold = $(x.thresholds[i])):")
        # println(io, format_confusion_matrix(x.confusion_matrices[i]))
    end
end

"""
    PolySurvivalFunction{T<:AbstractFloat}

A callable struct that represents a coral survival model using a regression.

# Fields
- `min_x::T` : Minimum diameter in training data
- `min_y::T` : Growth at minimum diameter
- `max_x::T` : Maximum diameter in training data
- `max_y::T` : Growth at maximum diameter
- `poly::Polynomial` : Polynomial
"""
struct PolySurvivalFunction{T<:AbstractFloat} <: Function
    min_x::T
    min_y::T
    max_x::T
    max_y::T
    poly::Polynomial

    function PolySurvivalFunction(xi::Vector{T}, yi::Vector{T}, poly::Polynomial) where T<:AbstractFloat
        new{T}(xi[1], yi[1], xi[end], yi[end], poly)
    end
end

# Make it callable using complementary log-log link
function (f::PolySurvivalFunction)(x::T) where T<:AbstractFloat
    return clamp(f.poly(log(x)), 0.0, 1.0)
end


"""
    PolySurvivalModel <: AbstractCoralBehavior

Functional relationships between coral size and survival.
"""
struct PolySurvivalModel <: AbstractCoralBehavior
    names::Vector{String}
    models::Vector{Function}
    performance::NamedTuple
end

function Base.show(io::IO, ::MIME"text/plain", x::PolySurvivalModel)
    title = "\nSurvival Model Performance Metrics"
    println(io, title)
    println(io, "─"^length(title))

    pretty_table(
        io,
        hcat(x.names, getfield.(x.models, :poly));
        header=["Group", "Model"]
    )

    explainer = """\n
        RMSE: 0.0 - Inf; Lower is better.
        R²: -∞ - 1.0; Higher is better.
        """
    println(io, explainer)

    data = hcat(
        x.names,
        x.performance.train.rmse,
        x.performance.test.rmse,
        x.performance.train.r2,
        x.performance.test.r2,
    )
    pretty_table(
        io, data;
        header=["Group", "Train RMSE", "Test RMSE", "Train R²", "Test R²"]
    )
end

# Load default survival models
survival_models = deserialize(joinpath(ASSET_DIR, "models", "survival_models.dat"))
