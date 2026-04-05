const _pi_f32 = Float32(π)
const _cover_scale = 1.0f0 / 10000.0f0

"""
    susceptibility_size_thresholds()::Vector{Float32}

Returns the minimum diameter (cm) at which corals in each functional group become 
susceptible to bleaching.
"""
function susceptibility_size_thresholds()::Vector{Float32}
    return Float32[
        12.5f0,  # Tabular Acropora
        12.5f0,  # Corymbose Acropora
        12.5f0,  # branching non-Acropora
        15.0f0,  # Small massives and encrusting
        20.0f0   # Large massives
    ]
end

"""
    mature_size_thresholds()::Vector{Float32}

Returns the minimum diameter (cm) at which corals in each functional group become 
reproductive.
"""
function mature_size_thresholds()::Vector{Float32}
    return Float32[
        5.0f0,  # Tabular Acropora
        5.0f0,  # Corymbose Acropora
        5.0f0,  # branching non-Acropora
        5.0f0,  # Small massives and encrusting
        5.0f0   # Large massives
    ]
end

"""
    bin_edges()::Matrix{Float32}

Helper function defining coral colony diameter bin edges in centimeters.
"""
function bin_edges()::Matrix{Float32}
    return Matrix{Float32}(
        [
            2.5f0 7.5f0 12.5f0 25.0f0 50.0f0 80.0f0 120.0f0 160.0f0;
            2.5f0 7.5f0 12.5f0 20.0f0 30.0f0 60.0f0 100.0f0 150.0f0;
            2.5f0 7.5f0 12.5f0 20.0f0 30.0f0 40.0f0 50.0f0 60.0f0;
            2.5f0 5.0f0 7.5f0 10.0f0 20.0f0 40.0f0 50.0f0 100.0f0;
            2.5f0 5.0f0 7.5f0 10.0f0 20.0f0 40.0f0 50.0f0 100.0f0
        ]
    )
end

"""
    diameter_size_classes()::Vector{Matrix{Float32}}

Determine diameter widths for each size class.

See also:
- `bin_edges()`
- `bin_widths()`
"""
function diameter_size_classes()::Vector{Matrix{Float32}}
    edges = bin_edges()
    n_groups = size(edges, 1)
    edges = hcat(zeros(Float32, n_groups), edges)

    # TODO: Use YAXArrays (groups ⋅ class start ⋅ class end)
    return map(x -> [x[1:(end - 1)] x[2:end]], eachrow(edges))
end

"""
    bin_widths()

Helper function defining coral colony diameter bin widths.
"""
function bin_widths()
    return bin_edges()[:, 2:end] .- bin_edges()[:, 1:(end - 1)]
end

function class_area()
    edges = bin_edges()
    return cover_cm_to_m2(edges)
end

"""
    linear_extensions()::Matrix{Float32}

Linear extensions. Data is the mean of functional group derived from ecoRRAP observations.
"""
function linear_extensions()::Matrix{Float32}
    return [
        2.45124f0 5.07098f0 5.01524f0 6.37291f0 6.72375f0 7.79938f0 0.0f0;
        2.4296f0 2.86086f0 2.78468f0 2.85766f0 3.14185f0 3.64447f0 0.0f0;
        1.75738f0 1.68217f0 1.50495f0 1.65701f0 1.65701f0 1.65701f0 0.0f0;
        1.16048f0 0.747208f0 0.748131f0 0.942616f0 1.33995f0 1.41176f0 0.0f0;
        1.1934f0 0.747208f0 0.748131f0 0.942616f0 1.33995f0 1.41176f0 0.0f0
    ]
end

"""
    survival_rates()::Matrix{Float32}

Survival rates. Data is mean of functional group derived from ecoRRAP observations.
"""
function survival_rates()::Matrix{Float32}
    return [
        0.687339f0 0.805556f0 0.788961f0 0.807143f0 0.842105f0 0.857143f0 0.857143f0  # Tabular Acropora
        0.776153f0 0.869252f0 0.908462f0 0.876652f0 0.889706f0 0.889706f0 0.889706f0  # Corymbose Acropora
        0.781176f0 0.871429f0 0.921466f0 0.916667f0 0.916667f0 0.916667f0 0.916667f0  # branching non-Acropora
        0.761658f0 0.920049f0 0.955396f0 0.973613f0 0.986486f0 0.984f0 0.972789f0  # Small massives and encrusting
        0.717391f0 0.920049f0 0.955396f0 0.973613f0 0.986486f0 0.984f0 0.972789f0  # Large massives
    ]
end

"""
Obtain the sum of all area (in m^2) for a given set of diameters.
"""
function cover_cm_to_m2(diameters::AbstractVector{F})::F where {F<:Float32}
    result = 0.0f0

    # Use @simd for vectorization, with a fused operation to reduce memory access
    @simd for d in diameters
        result += cover_cm_to_m2(d)
    end

    return result
end

"""
Convert centimeter diameter to meters area (m^2)
"""
function cover_cm_to_m2(d::F)::F where {F<:Float32}
    return _pi_f32 * (d * d * 0.25f0) * _cover_scale
end
function cover_cm_to_m2(
    diameters::AbstractArray{T}
)::AbstractArray{T} where {T<:AbstractFloat}
    return cover_cm_to_m2.(diameters)
end
function cover_cm_to_m2!(
    diameters::AbstractArray{T}, cache::AbstractArray{T}
)::Nothing where {T<:AbstractFloat}
    cache .= cover_cm_to_m2.(diameters)

    return nothing
end

"""
Convert area in m^2 to cm^2
"""
function cover_m2_to_cm(area::AbstractFloat)::AbstractFloat
    return 2.0f0 * sqrt(area / _pi_f32) * 100.0f0
end
