"""
    breeders(μ_o::T, μ_s::T, h²::T)::T where {T<:Float64}

Apply Breeder's equation.

```
S = μ_s - μ_o
μ_t1 = μ_o + (S * h²)
```

# Arguments
- `μ_o` : Mean of original population
- `μ_s` : Mean of next generation
- `h²` : Narrow-sense heritability

# Returns
Updated distribution mean
"""
function breeders(μ_o::T, μ_s::T, h²::T)::T where {T<:AbstractFloat}
    return μ_o + ((μ_s - μ_o) * h²)
end

"""
    larval_production(diam::Float32, grp::Int)::Float32

Calculate number of larvae produced by a coral colony based on its diameter.

# Extended help
Uses a power law relationship: `larvae = a * diameter^b`
where a and b are group-specific parameters.

Parameters tuned to produce reasonable recruitment rates based on colony size, but are
assumed values.

The method considers whole of life cycle, so the returned values are much lower than what
might be expected for egg production.
"""
function larval_production(diam::Float32, grp::Int)::Float32
    # Basic power law parameters for each group
    # a controls base fecundity, b controls how quickly it increases with size
    a = Float32[
        0.5f0,   # Tabular Acropora
        0.4f0,   # Corymbose Acropora
        0.3f0,   # Corymbose non-Acropora
        0.2f0,   # Small massives
        0.15f0   # Large massives
    ]

    # Hall and Hughes (1996) found production roughly scales with surface area
    b = Float32[2.0f0, 2.0f0, 2.0f0, 2.0f0, 2.0f0]

    return a[grp] * diam^b[grp]
end
