"""
    RMSE(y_hat, y)::Float64

Compute the Root Mean Square Error between model predictions and observed values.

RMSE measures the typical magnitude of prediction error in the same units as `y`.
For coral growth models predicting annual linear extension, an RMSE below 0.2 cm/yr
is generally considered a good fit given natural variability in EcoRRAP survey data.
For survival probability models (values in [0, 1]), RMSE below 0.1 is desirable.

\$\$
\\text{RMSE} = \\sqrt{\\frac{1}{n}\\sum_{i=1}^{n}(\\hat{y}_i - y_i)^2}
\$\$

# Arguments
- `y_hat` : Vector of model-predicted values (e.g., predicted annual linear extension
            in cm, or predicted survival probability).
- `y`     : Vector of observed field-survey values from EcoRRAP data, same length
            as `y_hat`.

# Returns
A non-negative `Float64` in the same units as `y`. A value of `0.0` indicates a
perfect fit.

# Notes
- **Units:** RMSE is expressed in the same units as `y`, making it directly
  interpretable (e.g., cm/yr for growth, dimensionless for survival probability).
- **Outlier sensitivity:** Squaring the residuals gives disproportionate weight to
  large errors; a single anomalous survey record can substantially inflate RMSE.
  Consider pairing with [`spearman`](@ref) or [`kendall`](@ref) for a
  rank-based sanity check.
- **Scale dependence:** RMSE should not be compared across datasets with different
  response scales (e.g., do not compare growth-model RMSE with survival-model RMSE).

# Examples
```jldoctest
julia> using Kora

julia> RMSE([0.5, 1.0, 1.5], [0.6, 0.9, 1.4])
0.1

julia> RMSE([1.0, 1.0, 1.0], [1.0, 1.0, 1.0])
0.0
```

# See Also
[`R2`](@ref), [`pearson`](@ref)
"""
RMSE(y_hat, y) = sqrt(mean((y_hat .- y) .^ 2))

"""
    R2(y_hat, y)::Float64

Compute the coefficient of determination (R-squared) between model predictions and
observed values.

R-squared quantifies the proportion of variance in the observed data that is explained
by the model. In EcoRRAP model evaluation, values above 0.8 indicate a strong fit;
values below 0.5 suggest the polynomial model captures less than half of the observed
variability, warranting model revision or additional covariates.

\$\$
R^2 = 1 - \\frac{SS_{res}}{SS_{tot}}, \\quad
SS_{res} = \\sum_{i}(y_i - \\hat{y}_i)^2, \\quad
SS_{tot} = \\sum_{i}(y_i - \\bar{y})^2
\$\$

# Arguments
- `y_hat` : Vector of model-predicted values (e.g., predicted annual linear extension
            in cm, or predicted survival probability).
- `y`     : Vector of observed field-survey values from EcoRRAP data, same length
            as `y_hat`.

# Returns
A `Float64`. A value of `1.0` indicates a perfect fit; `0.0` means the model
performs no better than predicting the mean of `y`; negative values indicate the
model performs worse than a constant mean predictor.

# Notes
- **Negative values:** R^2 is unbounded below zero. A negative R^2 signals a
  seriously misspecified model -- the polynomial may be predicting in the wrong
  direction or the training and evaluation datasets have incompatible distributions.
- **Undefined variance:** R^2 is undefined (division by zero) when all observed
  values are identical (i.e., `var(y) = 0`). This can occur with synthetic or
  heavily binned survey data.
- **Cross-dataset comparison:** R^2 values should not be compared between growth
  and survival models, or between different species groups, because SS_tot differs
  across response variables with different natural variances.

# Examples
```jldoctest
julia> using Kora

julia> R2([1.0, 2.0, 3.0], [1.0, 2.0, 3.0])
1.0

julia> R2([2.0, 2.0, 2.0], [1.0, 2.0, 3.0])
0.0
```

# See Also
[`RMSE`](@ref), [`pearson`](@ref)
"""
R2(y_hat, y) = 1 - (sum((y .- y_hat) .^ 2) / sum((y .- mean(y)) .^ 2))

"""
    pearson(y_hat, y)::Float64

Compute the Pearson product-moment correlation coefficient between model predictions
and observed values.

The Pearson correlation measures the strength and direction of the **linear** association
between `y_hat` and `y`. In coral reef model evaluation, a Pearson r above 0.9
indicates that the polynomial captures the dominant linear trend in the data, while
values below 0.7 suggest a weak linear relationship and possible model misspecification.

\$\$
r = \\frac{\\text{cov}(\\hat{y},\\, y)}{\\text{std}(\\hat{y}) \\cdot \\text{std}(y)}
\$\$

Implemented via `StatsBase.cor`.

# Arguments
- `y_hat` : Vector of model-predicted values (e.g., predicted annual linear extension
            in cm, or predicted survival probability).
- `y`     : Vector of observed field-survey values from EcoRRAP data, same length
            as `y_hat`.

# Returns
A `Float64` in [-1, 1]. A value of `1.0` indicates a perfect positive linear
relationship; `-1.0` indicates a perfect negative linear relationship; `0.0`
indicates no linear association.

# Notes
- **Linearity assumption:** Pearson r only captures **linear** relationships. Coral
  growth responses are often nonlinear across size classes; a low Pearson r alongside
  a high [`spearman`](@ref) r would indicate a monotonic but nonlinear relationship.
- **Outlier sensitivity:** Like [`RMSE`](@ref), Pearson r is sensitive to extreme
  survey records. A single anomalous data point can substantially shift the estimate.
- **Scale invariance:** Unlike [`RMSE`](@ref), Pearson r is dimensionless and
  comparable across growth and survival models.

# Examples
```jldoctest
julia> using Kora

julia> pearson([1.0, 2.0, 3.0], [1.0, 2.0, 3.0])
1.0

julia> pearson([1.0, 2.0, 3.0], [-1.0, -2.0, -3.0])
-1.0
```

# See Also
[`spearman`](@ref), [`kendall`](@ref), [`R2`](@ref)
"""
pearson(y_hat, y) = StatsBase.cor(y_hat, y)

"""
    spearman(y_hat, y)::Float64

Compute the Spearman rank correlation coefficient between model predictions and
observed values.

Spearman's rho measures the strength of the **monotonic** (not necessarily linear)
relationship between `y_hat` and `y` by operating on the ranks of the data rather
than the raw values. In EcoRRAP model evaluation it is preferable to [`pearson`](@ref)
when the growth or survival response is nonlinear, or when the survey dataset contains
outlier records that cannot be cleaned before evaluation.

\$\$
\\rho = r\\bigl(\\text{rank}(\\hat{y}),\\, \\text{rank}(y)\\bigr)
\$\$

Implemented via `StatsBase.corspearman`.

# Arguments
- `y_hat` : Vector of model-predicted values (e.g., predicted annual linear extension
            in cm, or predicted survival probability).
- `y`     : Vector of observed field-survey values from EcoRRAP data, same length
            as `y_hat`.

# Returns
A `Float64` in [-1, 1]. A value of `1.0` indicates a perfect positive monotonic
relationship; `-1.0` indicates a perfect negative monotonic relationship; `0.0`
indicates no monotonic association.

# Notes
- **Outlier robustness:** Because ranks are bounded by [1, n], extreme survey
  values inflate Spearman rho far less than they inflate [`pearson`](@ref) or
  [`RMSE`](@ref).
- **Nonlinear sensitivity:** Spearman rho captures any monotonic relationship, so
  a model that preserves rank order (even with a nonlinear bias) will score well.
  Pair with [`RMSE`](@ref) to distinguish rank-correct-but-biased models.
- **Tied ranks:** Ties in `y` or `y_hat` are resolved using the average-rank
  convention in `StatsBase.corspearman`.

# Examples
```jldoctest
julia> using Kora

julia> spearman([1.0, 2.0, 3.0], [4.0, 5.0, 6.0])
1.0

julia> spearman([1.0, 2.0, 3.0], [6.0, 5.0, 4.0])
-1.0
```

# See Also
[`pearson`](@ref), [`kendall`](@ref)
"""
spearman(y_hat, y) = StatsBase.corspearman(y_hat, y)

"""
    kendall(y_hat, y)::Float64

Compute Kendall's tau-b rank correlation coefficient between model predictions and
observed values.

Kendall's tau measures **rank concordance**: the proportion of observation pairs
whose relative ordering is consistent between `y_hat` and `y`, minus the proportion
that is discordant. It is more conservative than [`spearman`](@ref) and is preferred
for small EcoRRAP subgroup evaluations (e.g., per-species or per-reef-zone subsets)
because its sampling distribution is better-behaved at small n.

\$\$
\\tau_b = \\frac{C - D}{\\sqrt{(C + D + T_x)(C + D + T_y)}}
\$\$

where ``C`` is the number of concordant pairs, ``D`` the number of discordant pairs,
and ``T_x``, ``T_y`` are tie counts in `y_hat` and `y` respectively.
Implemented via `StatsBase.corkendall`.

# Arguments
- `y_hat` : Vector of model-predicted values (e.g., predicted annual linear extension
            in cm, or predicted survival probability).
- `y`     : Vector of observed field-survey values from EcoRRAP data, same length
            as `y_hat`.

# Returns
A `Float64` in [-1, 1]. A value of `1.0` indicates perfect concordance (all pairs
in the same order); `-1.0` indicates perfect discordance; `0.0` indicates no
rank association.

# Notes
- **Small-sample robustness:** Kendall's tau has a more tractable exact null
  distribution than [`spearman`](@ref) when n < 30, which is common in per-reef
  or per-species subgroup evaluations.
- **Computational cost:** Kendall's tau requires O(n^2) pair comparisons versus
  O(n log n) for Spearman. For the full EcoRRAP dataset this difference is
  negligible, but it may matter if `kendall` is called inside a large bootstrap
  loop.
- **Interpretation vs. Spearman:** |tau| is numerically smaller than |rho| for
  the same data; the two metrics should not be compared by magnitude across
  functions -- use [`ALL_METRICS`](@ref) to report them side-by-side.

# Examples
```jldoctest
julia> using Kora

julia> kendall([1.0, 2.0, 3.0], [1.0, 2.0, 3.0])
1.0

julia> kendall([1.0, 2.0, 3.0], [3.0, 2.0, 1.0])
-1.0
```

# See Also
[`spearman`](@ref), [`pearson`](@ref)
"""
kendall(y_hat, y) = StatsBase.corkendall(y_hat, y)

"""
    ALL_METRICS

Ordered vector of all five performance metric functions used to evaluate fitted
polynomial growth and survival models against EcoRRAP held-out field data.

Contains, in order: [`RMSE`](@ref), [`R2`](@ref), [`pearson`](@ref),
[`spearman`](@ref), [`kendall`](@ref).

Each element is a callable with the shared signature `f(y_hat, y)`, where `y_hat`
is the vector of model-predicted values and `y` is the vector of observed values.
This uniform signature allows `ALL_METRICS` to be iterated directly when populating
the performance table in `process_ecorrap_models`, and ensures the table columns
produced by `show(::PolyGrowthModel)` and `show(::PolySurvivalModel)` are always
in a consistent order.

# Examples
```jldoctest
julia> using Kora

julia> length(ALL_METRICS)
5

julia> all(f([1.0, 2.0, 3.0], [1.0, 2.0, 3.0]) isa Float64 for f in ALL_METRICS)
true
```

# See Also
[`RMSE`](@ref), [`R2`](@ref), [`pearson`](@ref), [`spearman`](@ref),
[`kendall`](@ref)
"""
const ALL_METRICS = [RMSE, R2, pearson, spearman, kendall]
