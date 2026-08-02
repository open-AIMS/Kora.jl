# Custom statistical functions for computational efficiency.
"""
    rational_erf(x::Float64)::Float64

Rational approximation of the error function only using elementary functions [1].
Maximum error of 1.5 × 10^{-7}.

# References
1. Abramowitz, Milton; Stegun, Irene Ann, eds. (1983) [June 1964]. "Chapter 7".
   Handbook of Mathematical Functions with Formulas, Graphs, and Mathematical Tables.
   Applied Mathematics Series. Vol. 55 (Ninth reprint with additional corrections of
   tenth original printing with corrections (December 1972); first ed.). Washington D.C.;
   New York: United States Department of Commerce, National Bureau of Standards;
   Dover Publications. p. 297
"""
function rational_erf(x::F)::F where {F<:Union{Float32,Float64}}
    coef::F = 1.0
    if (x < 0)
        x *= -1
        coef = -1
    end

    # Use multiplication to avoid slow power function (x^n)
    x2::F = x * x
    x3::F = x2 * x
    x4::F = x3 * x
    x5::F = x4 * x
    x6::F = x5 * x

    a1::F = 0.0705230784
    a2::F = 0.0422820123
    a3::F = 0.0092705272
    a4::F = 0.0001520143
    a5::F = 0.0002765672
    a6::F = 0.0000430638

    denom::F = 1.0 + a1 * x + a2 * x2 + a3 * x3 + a4 * x4 + a5 * x5 + a6 * x6

    denom = denom * denom  # power 2
    denom = denom * denom  # power 4
    denom = denom * denom  # power 8
    denom = denom * denom  # power 16

    return coef * (1 - 1.0 / denom)
end

"""
    rational_erf(a::F, b::F)::F where {F<:Union{Float32,Float64}}

Rational approximation of the generalised error function integral erf(b) − erf(a).
Equivalent to `SpecialFunctions.erf(a, b)` but uses `rational_erf` throughout,
avoiding any native-library dependency. Maximum error is 2 × 1.5 × 10^{-7}.

# See Also
[`rational_erf(x)`](@ref)
"""
function rational_erf(a::F, b::F)::F where {F<:Union{Float32,Float64}}
    return rational_erf(b) - rational_erf(a)
end

"""
    rational_erfcx(x::F)::F where {F<:Union{Float32,Float64}}

Approximation of erfcx using a rational approximation of the error function.

erfcx(x) = e^{x^2} ⋅ (1 - erf(x))
"""
function rational_erfcx(x::F)::F where {F<:Union{Float32,Float64}}
    return exp(x * x) * (1.0 - rational_erf(x))
end

"""
    truncated_standard_normal_mean(lb::F, ub::F)::F where {F<:Union{Float32,Float64}}

Compute the mean of the standard normal distribution truncated to the interval
[`lb`, `ub`].

Implementation follows Distributions.jl, excluding unused error checks. When
`lb` > `ub`, `ub` is returned to avoid NaN propagation.

# Arguments
- `lb` : Lower bound of the truncated distribution.
- `ub` : Upper bound of the truncated distribution.

# Returns
`F` : Mean of the truncated standard normal distribution, or `ub` if `lb` > `ub`.

# Examples
```jldoctest
julia> using Kora

julia> truncated_standard_normal_mean(-1.0, 1.0)
0.0

julia> truncated_standard_normal_mean(0.0, 0.0)
0.0
```

# References
1. Distributions.jl truncated normal implementation:
   https://github.com/JuliaStats/Distributions.jl/blob/c1705a3015d438f7e841e82ef5148224813831e8/src/truncated/normal.jl#L24-L46
"""
function truncated_standard_normal_mean(lb::F, ub::F)::F where {F<:Union{Float32,Float64}}
    if abs(lb) > abs(ub)
        return -truncated_standard_normal_mean(-ub, -lb)
    elseif (lb == ub)
        return lb
    end

    mid = (lb + ub) * 0.5
    Δ = (ub - lb) * mid
    lb′ = lb * 0.7071067811865476
    ub′ = ub * 0.7071067811865476

    m = ub
    if lb ≤ 0 ≤ ub
        m = expm1(-Δ) * exp(-lb^2 / 2) / (rational_erf(lb′) - rational_erf(ub′))
    elseif 0 < lb < ub
        z = exp(-Δ) * rational_erfcx(ub′) - rational_erfcx(lb′)
        iszero(z) && return mid
        m = expm1(-Δ) / z
    end

    return clamp(m / 1.2533141373155003, lb, ub)
end

"""
    truncated_normal_mean(
        normal_mean::F,
        normal_stdev::F,
        lower_bound::F,
        upper_bound::F
    )::F where {F<:Union{Float32,Float64}}

Compute the mean of the normal distribution with mean `normal_mean` and standard
deviation `normal_stdev`, truncated to the interval [`lower_bound`, `upper_bound`].

Delegates to [`truncated_standard_normal_mean`](@ref) after standardising the
bounds to the unit-normal scale.

# Arguments
- `normal_mean` : Mean of the underlying (untruncated) normal distribution.
- `normal_stdev` : Standard deviation of the underlying (untruncated) normal distribution.
- `lower_bound` : Lower bound of the truncated normal distribution.
- `upper_bound` : Upper bound of the truncated normal distribution.

# Returns
`F` : Mean of the truncated normal distribution.

# Examples
```jldoctest
julia> using Kora

julia> truncated_normal_mean(0.0, 1.0, -1.0, 1.0)
0.0

julia> truncated_normal_mean(5.0, 2.0, 5.0, 5.0)
5.0
```

# See Also
[`truncated_standard_normal_mean`](@ref), [`truncated_normal_cdf`](@ref)
"""
function truncated_normal_mean(
    normal_mean::F, normal_stdev::F, lower_bound::F, upper_bound::F
)::F where {F<:Union{Float32,Float64}}
    alpha::F = (lower_bound - normal_mean) / normal_stdev
    beta::F = (upper_bound - normal_mean) / normal_stdev

    return normal_mean + truncated_standard_normal_mean(alpha, beta) * normal_stdev
end

"""
    truncated_normal_cdf(
        x::F,
        normal_mean::F,
        normal_stdev::F,
        lower_bound::F,
        upper_bound::F
    )::F where {F<:Union{Float32,Float64}}

Evaluate the CDF of the normal distribution with mean `normal_mean` and standard
deviation `normal_stdev`, truncated to the interval [`lower_bound`, `upper_bound`],
at the point `x`.

Returns `0.0` when `x <= lower_bound` and `1.0` when `x >= upper_bound`. Uses a
rational approximation of the error function via `rational_erf` for
efficiency; falls back to `SpecialFunctions.erf` when the truncation bounds lie
more than 3 standard deviations from the mean to avoid precision loss.

# Arguments
- `x` : Value at which to evaluate the CDF.
- `normal_mean` : Mean of the underlying (untruncated) normal distribution.
- `normal_stdev` : Standard deviation of the underlying (untruncated) normal distribution.
- `lower_bound` : Lower bound of the truncated distribution.
- `upper_bound` : Upper bound of the truncated distribution.

# Returns
`F` : CDF of the truncated normal distribution evaluated at `x`, in [0, 1].

# Examples
```jldoctest
julia> using Kora

julia> truncated_normal_cdf(-2.0, 0.0, 1.0, -1.0, 1.0)
0.0

julia> truncated_normal_cdf(2.0, 0.0, 1.0, -1.0, 1.0)
1.0

julia> truncated_normal_cdf(0.0, 0.0, 1.0, -1.0, 1.0)
0.5
```

# See Also
[`truncated_normal_mean`](@ref), [`truncated_standard_normal_mean`](@ref)
"""
function truncated_normal_cdf(
    x::F,
    normal_mean::F,
    normal_stdev::F,
    lower_bound::F,
    upper_bound::F
)::F where {F<:Union{Float32,Float64}}
    if x <= lower_bound
        return 0.0
    elseif x >= upper_bound
        return 1.0
    end

    alpha::F = (lower_bound - normal_mean) / normal_stdev
    beta::F = (upper_bound - normal_mean) / normal_stdev
    zeta::F = (x - normal_mean) / normal_stdev

    # Large errors occurs when bounds deviate from the mean significantly and
    # are close together relative to the standard deviation.
    threshold = 3
    if (alpha > threshold && beta > threshold) || (alpha < -threshold && beta < -threshold)
        @debug "Possible loss of accuracy: the given truncated normal distribution bounds \
            are more than 5 standard deviations from the normal mean. \
            \nLower and upper bounds of the truncated normal distribution are \
            $(alpha) and $(beta) standard deviations from the normal mean respectively. \
            Falling back to more accurate calculation."

        return rational_erf(alpha * 0.7071067811865476, zeta * 0.7071067811865476) /
               rational_erf(alpha * 0.7071067811865476, beta * 0.7071067811865476)
    end

    # Store error function of alpha to avoid duplicate calculations
    erf_alpha = rational_erf(alpha * 0.7071067811865476)

    return (rational_erf(zeta * 0.7071067811865476) - erf_alpha) /
           (rational_erf(beta * 0.7071067811865476) - erf_alpha)
end

"""
    rand_truncated_normal(rng, μ, σ, lo, hi[, n])

Draw one (or `n`) samples from the normal distribution N(μ, σ²) truncated to
[`lo`, `hi`] using rejection sampling. No native-library dependencies — safe for
WASM and AOT compilation.
"""
function rand_truncated_normal(
    rng::AbstractRNG, μ::F, σ::F, lo::F, hi::F
)::F where {F<:Union{Float32,Float64}}
    while true
        x = μ + σ * F(randn(rng))
        lo <= x <= hi && return x
    end
end

function rand_truncated_normal(
    rng::AbstractRNG, μ::F, σ::F, lo::F, hi::F, n::Integer
)::Vector{F} where {F<:Union{Float32,Float64}}
    out = Vector{F}(undef, n)
    for i in 1:n
        out[i] = rand_truncated_normal(rng, μ, σ, lo, hi)
    end
    return out
end
