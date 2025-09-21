# Private consts - only used in this file
const _euler_f32b = Float32(ℯ)

"""
    growth_inflection_point()::Vector{Float32}

Value indicates when to begin adjusting growth rates as cover approaches the maximum
available hard substrate.

See: `space_constraint()`
"""
function growth_inflection_point()::Vector{Float32}
    return Float32[0.96, 0.96, 0.96, 0.96, 0.96, 0.96]
end

"""
    space_constraint(x::F, k::F; x0::F=0.96f0)::F where {F<:Float32}

Modify growth as coral cover approaches total habitable area.

# Arguments
- `x` : proportion of available space currently taken
- `k` : steepness parameter, where the larger the value, the faster the drop (suggest k=20).
- `x0` : Inflection point where growth begins to be constrained (default: 96% of available area).
"""
@inline function space_constraint(x::F, k::F; x0::F=0.96f0)::F where {F<:Float32}
    if x >= 1.0f0
        return 0.0f0
    end

    return 1.0f0 / (1.0f0 + _euler_f32b^(k * (x - x0)))
end

# Growth models need to account for available space... (see space constraint)
function growth!(
    model::M,
    diam::Vector{Vector{F}},
    reef_cover::Vector{F},
    grp_mod::F,
    loc_scalers::Vector{F}
)::Nothing where {M,F<:Float32}
    # Use explicit loops to avoid broadcast allocations
    Threads.@threads for i in eachindex(reef_cover)
        constraint = space_constraint(reef_cover[i], 20.0f0; x0=grp_mod)
        scaler = loc_scalers[i]

        for j in eachindex(diam[i])
            old_diam = diam[i][j]
            diam[i][j] = old_diam + (model(old_diam) - old_diam) * constraint * scaler
        end
    end

    return nothing
end

struct PolyGrowthModel <: AbstractCoralBehavior
    names::Vector{String}
    models::Vector{Function}
    performance::NamedTuple
end

function Base.show(io::IO, ::MIME"text/plain", x::PolyGrowthModel)
    title = "\nGrowth Model Performance Metrics"
    println(io, title)
    println(io, "─"^length(title))

    pretty_table(
        io,
        hcat(x.names, getfield.(x.models, :poly));
        header=["Group", "Model"]
    )

    explainer = """
        RMSE: 0.0 - Inf; Lower is better.
        R²: -∞ - 1.0; Higher is better.
        """
    println(io, explainer)

    performance_matrix = [
        (getfield(x.performance.train, m), getfield(x.performance.test, m))
        for m in Symbol.(CoralFlow.ALL_METRICS)
    ]

    performances = hcat([hcat(t[1], t[2]) for t in performance_matrix]...)
    perf_heads = [("Train $m", "Test $m") for m in Symbol.(ALL_METRICS)]
    perf_headers = vcat([vcat(t[1], t[2]) for t in perf_heads]...)

    data = hcat(
        x.names,
        performances
    )
    pretty_table(
        io, data;
        header=["Group", perf_headers...]
    )
end

"""
    PolyGrowthFunction{T<:AbstractFloat}

A callable struct that represents a coral growth model.

# Fields
- `min_x::T` : Minimum diameter in training data
- `min_y::T` : Growth at minimum diameter
- `max_x::T` : Maximum diameter in training data
- `max_y::T` : Growth at maximum diameter
- `poly::Polynomial` : Polynomial
"""
struct PolyGrowthFunction{T<:AbstractFloat} <: Function
    min_x::T
    min_y::T
    max_x::T
    max_y::T
    poly::Polynomial

    function PolyGrowthFunction(xi::Vector{T}, yi::Vector{T}, poly::Polynomial) where T<:AbstractFloat
        new{T}(xi[1], yi[1], xi[end], yi[end], poly)
    end
end

"""Define `PolyGrowthFunction` as a callable."""
function (f::PolyGrowthFunction)(x::T) where T<:AbstractFloat
    if x < f.min_x
        return f.min_y
    elseif x > f.max_x
        return f.max_y
    end

    y = exp(f.poly(log(x)))
    if y > f.max_y
        return f.max_y
    end

    return y
end

# Load default growth models
growth_models = deserialize(joinpath(ASSET_DIR, "models", "Offshore_Northern_growth_models.dat"))
