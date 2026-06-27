import DataFrames.PrettyTables: pretty_table
using Polynomials: Polynomial

# Private consts - only used in this file
const _euler_f32b = Float32(ℯ)

"""
    growth_inflection_point()::Vector{Float32}

Value indicates when to begin adjusting growth rates as cover approaches the maximum
available hard substrate.

See: `space_constraint()`
"""
function growth_inflection_point()::Vector{Float32}
    return Float32[0.96f0, 0.96f0, 0.96f0, 0.96f0, 0.96f0, 0.96f0]
end

"""
    space_constraint(x::F, k::F; x0::F=0.96f0)::F where {F<:Float32}

Modify growth as coral cover approaches total habitable area.

# Arguments
- `x` : proportion of available space currently taken
- `k` : steepness parameter, where the larger the value, the faster the drop (recommended: k=20).
- `x0` : Inflection point where growth begins to be constrained (default: 96% of available area).
"""
@inline function space_constraint(x::F, k::F; x0::F=0.96f0)::F where {F<:Float32}
    if x >= 1.0f0
        return 0.0f0
    end

    return 1.0f0 / (1.0f0 + _euler_f32b^(k * (x - x0)))
end

"""
    PolyGrowthFunction{T<:AbstractFloat,P<:Polynomial} <: Function

A callable struct that represents a coral growth model.

# Fields
- `min_x::T` : Minimum diameter in training data
- `min_y::T` : Growth at minimum diameter
- `max_x::T` : Maximum diameter in training data
- `max_y::T` : Growth at maximum diameter
- `poly::P` : Polynomial
"""
struct PolyGrowthFunction{T<:AbstractFloat,P<:Polynomial} <: Function
    min_x::T
    min_y::T
    max_x::T
    max_y::T
    poly::P

    function PolyGrowthFunction(
        xi::Vector{T}, yi::Vector{T}, poly::P
    ) where {T<:AbstractFloat,P<:Polynomial}
        return new{T,P}(xi[1], yi[1], xi[end], yi[end], poly)
    end

    function PolyGrowthFunction(
        xi::Vector{T},
        yi::Vector{T},
        poly::P,
        max_x::AbstractFloat,
        max_y::AbstractFloat
    ) where {T<:AbstractFloat,P<:Polynomial}
        return new{T,P}(
            Float32(xi[1]), Float32(yi[1]), Float32(max_x), Float32(max_y), poly
        )
    end

    # Bare-fields constructor for JSON round-trip loading.
    # All fields are provided directly; no xi/yi vector inference.
    # Must not reference any module-level mutable bindings.
    function PolyGrowthFunction(
        min_x::T, min_y::T, max_x::T, max_y::T, poly::P
    ) where {T<:AbstractFloat,P<:Polynomial}
        return new{T,P}(min_x, min_y, max_x, max_y, poly)
    end
end

"""Define `PolyGrowthFunction` as a callable."""
function (f::PolyGrowthFunction)(x::T)::T where T<:Float32
    if x < f.min_x
        return f.min_y
    end

    return min(f.poly(log(x)), f.max_y)
end

struct PolyGrowthModel <: AbstractCoralBehavior
    names::Vector{String}
    models::Vector{PolyGrowthFunction}
    performance::NamedTuple
end

Base.length(m::PolyGrowthModel) = length(m.models)
Base.getindex(m::PolyGrowthModel, i::Int) = m.models[i]

function Base.show(io::IO, ::MIME"text/plain", x::PolyGrowthModel)
    title = "\nGrowth Model Performance Metrics"
    println(io, title)
    println(io, "─"^length(title))

    pretty_table(
        io,
        hcat(x.names, getfield.(x.models, :poly));
        column_labels=["Group", "Model"]
    )

    explainer = """
        RMSE: 0.0 - Inf; Lower is better.
        R²: -∞ - 1.0; Higher is better.
        """
    println(io, explainer)

    performance_matrix = [
        (getfield(x.performance.train, m), getfield(x.performance.test, m))
        for m in Symbol.(Kora.ALL_METRICS)
    ]

    performances = hcat([hcat(t[1], t[2]) for t in performance_matrix]...)
    perf_heads = [("Train $m", "Test $m") for m in Symbol.(ALL_METRICS)]
    perf_headers = vcat([vcat(t[1], t[2]) for t in perf_heads]...)

    data = hcat(
        x.names,
        round.(performances; digits=3)
    )
    return pretty_table(
        io, data;
        column_labels=["Group", perf_headers...]
    )
end

growth_models::Union{Nothing,PolyGrowthModel} = nothing
