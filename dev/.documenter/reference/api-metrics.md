
# Metrics API {#Metrics-API}
<details class='jldocstring custom-block' open>
<summary><a id='Kora.RMSE' href='#Kora.RMSE'><span class="jlbinding">Kora.RMSE</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
RMSE(y_hat, y)::Float64
```


Compute the Root Mean Square Error between model predictions and observed values.

RMSE measures the typical magnitude of prediction error in the same units as `y`. For coral growth models predicting annual linear extension, an RMSE below 0.2 cm/yr is generally considered a good fit given natural variability in EcoRRAP survey data. For survival probability models (values in [0, 1]), RMSE below 0.1 is desirable.

$

\text{RMSE} = \sqrt{\frac{1}{n}\sum_{i=1}^{n}(\hat{y}_i - y_i)^2} $

**Arguments**
- `y_hat` : Vector of model-predicted values (e.g., predicted annual linear extension           in cm, or predicted survival probability).
  
- `y`     : Vector of observed field-survey values from EcoRRAP data, same length           as `y_hat`.
  

**Returns**

A non-negative `Float64` in the same units as `y`. A value of `0.0` indicates a perfect fit.

**Notes**
- **Units:** RMSE is expressed in the same units as `y`, making it directly interpretable (e.g., cm/yr for growth, dimensionless for survival probability).
  
- **Outlier sensitivity:** Squaring the residuals gives disproportionate weight to large errors; a single anomalous survey record can substantially inflate RMSE. Consider pairing with [`spearman`](/reference/api-metrics#Kora.spearman) or [`kendall`](/reference/api-metrics#Kora.kendall) for a rank-based sanity check.
  
- **Scale dependence:** RMSE should not be compared across datasets with different response scales (e.g., do not compare growth-model RMSE with survival-model RMSE).
  

**Examples**

```julia
julia> using Kora

julia> RMSE([0.5, 1.0, 1.5], [0.6, 0.9, 1.4])
0.1

julia> RMSE([1.0, 1.0, 1.0], [1.0, 1.0, 1.0])
0.0
```


**See Also**

[`R2`](/reference/api-metrics#Kora.R2), [`pearson`](/reference/api-metrics#Kora.pearson)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/Kora.jl/blob/9ccdfd81b778a9c3de77df920a4498af80a8832a/src/metrics.jl#L1-L48" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Kora.R2' href='#Kora.R2'><span class="jlbinding">Kora.R2</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
R2(y_hat, y)::Float64
```


Compute the coefficient of determination (R-squared) between model predictions and observed values.

R-squared quantifies the proportion of variance in the observed data that is explained by the model. In EcoRRAP model evaluation, values above 0.8 indicate a strong fit; values below 0.5 suggest the polynomial model captures less than half of the observed variability, warranting model revision or additional covariates.

$

R^2 = 1 - \frac{SS_{res}}{SS_{tot}}, \quad SS_{res} = \sum_{i}(y_i - \hat{y}_i)^2, \quad SS_{tot} = \sum_{i}(y_i - \bar{y})^2 $

**Arguments**
- `y_hat` : Vector of model-predicted values (e.g., predicted annual linear extension           in cm, or predicted survival probability).
  
- `y`     : Vector of observed field-survey values from EcoRRAP data, same length           as `y_hat`.
  

**Returns**

A `Float64`. A value of `1.0` indicates a perfect fit; `0.0` means the model performs no better than predicting the mean of `y`; negative values indicate the model performs worse than a constant mean predictor.

**Notes**
- **Negative values:** R^2 is unbounded below zero. A negative R^2 signals a seriously misspecified model – the polynomial may be predicting in the wrong direction or the training and evaluation datasets have incompatible distributions.
  
- **Undefined variance:** R^2 is undefined (division by zero) when all observed values are identical (i.e., `var(y) = 0`). This can occur with synthetic or heavily binned survey data.
  
- **Cross-dataset comparison:** R^2 values should not be compared between growth and survival models, or between different species groups, because SS_tot differs across response variables with different natural variances.
  

**Examples**

```julia
julia> using Kora

julia> R2([1.0, 2.0, 3.0], [1.0, 2.0, 3.0])
1.0

julia> R2([2.0, 2.0, 2.0], [1.0, 2.0, 3.0])
0.0
```


**See Also**

[`RMSE`](/reference/api-metrics#Kora.RMSE), [`pearson`](/reference/api-metrics#Kora.pearson)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/Kora.jl/blob/9ccdfd81b778a9c3de77df920a4498af80a8832a/src/metrics.jl#L51-L103" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Kora.pearson' href='#Kora.pearson'><span class="jlbinding">Kora.pearson</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
pearson(y_hat, y)::Float64
```


Compute the Pearson product-moment correlation coefficient between model predictions and observed values.

The Pearson correlation measures the strength and direction of the **linear** association between `y_hat` and `y`. In coral reef model evaluation, a Pearson r above 0.9 indicates that the polynomial captures the dominant linear trend in the data, while values below 0.7 suggest a weak linear relationship and possible model misspecification.

$

r = \frac{\text{cov}(\hat{y},\, y)}{\text{std}(\hat{y}) \cdot \text{std}(y)} $

Implemented via `StatsBase.cor`.

**Arguments**
- `y_hat` : Vector of model-predicted values (e.g., predicted annual linear extension           in cm, or predicted survival probability).
  
- `y`     : Vector of observed field-survey values from EcoRRAP data, same length           as `y_hat`.
  

**Returns**

A `Float64` in [-1, 1]. A value of `1.0` indicates a perfect positive linear relationship; `-1.0` indicates a perfect negative linear relationship; `0.0` indicates no linear association.

**Notes**
- **Linearity assumption:** Pearson r only captures **linear** relationships. Coral growth responses are often nonlinear across size classes; a low Pearson r alongside a high [`spearman`](/reference/api-metrics#Kora.spearman) r would indicate a monotonic but nonlinear relationship.
  
- **Outlier sensitivity:** Like [`RMSE`](/reference/api-metrics#Kora.RMSE), Pearson r is sensitive to extreme survey records. A single anomalous data point can substantially shift the estimate.
  
- **Scale invariance:** Unlike [`RMSE`](/reference/api-metrics#Kora.RMSE), Pearson r is dimensionless and comparable across growth and survival models.
  

**Examples**

```julia
julia> using Kora

julia> pearson([1.0, 2.0, 3.0], [1.0, 2.0, 3.0])
1.0

julia> pearson([1.0, 2.0, 3.0], [-1.0, -2.0, -3.0])
-1.0
```


**See Also**

[`spearman`](/reference/api-metrics#Kora.spearman), [`kendall`](/reference/api-metrics#Kora.kendall), [`R2`](/reference/api-metrics#Kora.R2)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/Kora.jl/blob/9ccdfd81b778a9c3de77df920a4498af80a8832a/src/metrics.jl#L106-L156" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Kora.spearman' href='#Kora.spearman'><span class="jlbinding">Kora.spearman</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
spearman(y_hat, y)::Float64
```


Compute the Spearman rank correlation coefficient between model predictions and observed values.

Spearman's rho measures the strength of the **monotonic** (not necessarily linear) relationship between `y_hat` and `y` by operating on the ranks of the data rather than the raw values. In EcoRRAP model evaluation it is preferable to [`pearson`](/reference/api-metrics#Kora.pearson) when the growth or survival response is nonlinear, or when the survey dataset contains outlier records that cannot be cleaned before evaluation.

$

\rho = r\bigl(\text{rank}(\hat{y}),\, \text{rank}(y)\bigr) $

Implemented via `StatsBase.corspearman`.

**Arguments**
- `y_hat` : Vector of model-predicted values (e.g., predicted annual linear extension           in cm, or predicted survival probability).
  
- `y`     : Vector of observed field-survey values from EcoRRAP data, same length           as `y_hat`.
  

**Returns**

A `Float64` in [-1, 1]. A value of `1.0` indicates a perfect positive monotonic relationship; `-1.0` indicates a perfect negative monotonic relationship; `0.0` indicates no monotonic association.

**Notes**
- **Outlier robustness:** Because ranks are bounded by [1, n], extreme survey values inflate Spearman rho far less than they inflate [`pearson`](/reference/api-metrics#Kora.pearson) or [`RMSE`](/reference/api-metrics#Kora.RMSE).
  
- **Nonlinear sensitivity:** Spearman rho captures any monotonic relationship, so a model that preserves rank order (even with a nonlinear bias) will score well. Pair with [`RMSE`](/reference/api-metrics#Kora.RMSE) to distinguish rank-correct-but-biased models.
  
- **Tied ranks:** Ties in `y` or `y_hat` are resolved using the average-rank convention in `StatsBase.corspearman`.
  

**Examples**

```julia
julia> using Kora

julia> spearman([1.0, 2.0, 3.0], [4.0, 5.0, 6.0])
1.0

julia> spearman([1.0, 2.0, 3.0], [6.0, 5.0, 4.0])
-1.0
```


**See Also**

[`pearson`](/reference/api-metrics#Kora.pearson), [`kendall`](/reference/api-metrics#Kora.kendall)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/Kora.jl/blob/9ccdfd81b778a9c3de77df920a4498af80a8832a/src/metrics.jl#L159-L211" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Kora.kendall' href='#Kora.kendall'><span class="jlbinding">Kora.kendall</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
kendall(y_hat, y)::Float64
```


Compute Kendall's tau-b rank correlation coefficient between model predictions and observed values.

Kendall's tau measures **rank concordance**: the proportion of observation pairs whose relative ordering is consistent between `y_hat` and `y`, minus the proportion that is discordant. It is more conservative than [`spearman`](/reference/api-metrics#Kora.spearman) and is preferred for small EcoRRAP subgroup evaluations (e.g., per-species or per-reef-zone subsets) because its sampling distribution is better-behaved at small n.

$

\tau_b = \frac{C - D}{\sqrt{(C + D + T_x)(C + D + T_y)}} $

where $C$ is the number of concordant pairs, $D$ the number of discordant pairs, and $T_x$, $T_y$ are tie counts in `y_hat` and `y` respectively. Implemented via `StatsBase.corkendall`.

**Arguments**
- `y_hat` : Vector of model-predicted values (e.g., predicted annual linear extension           in cm, or predicted survival probability).
  
- `y`     : Vector of observed field-survey values from EcoRRAP data, same length           as `y_hat`.
  

**Returns**

A `Float64` in [-1, 1]. A value of `1.0` indicates perfect concordance (all pairs in the same order); `-1.0` indicates perfect discordance; `0.0` indicates no rank association.

**Notes**
- **Small-sample robustness:** Kendall's tau has a more tractable exact null distribution than [`spearman`](/reference/api-metrics#Kora.spearman) when n &lt; 30, which is common in per-reef or per-species subgroup evaluations.
  
- **Computational cost:** Kendall's tau requires O(n^2) pair comparisons versus O(n log n) for Spearman. For the full EcoRRAP dataset this difference is negligible, but it may matter if `kendall` is called inside a large bootstrap loop.
  
- **Interpretation vs. Spearman:** |tau| is numerically smaller than |rho| for the same data; the two metrics should not be compared by magnitude across functions – use [`ALL_METRICS`](/reference/api-metrics#Kora.ALL_METRICS) to report them side-by-side.
  

**Examples**

```julia
julia> using Kora

julia> kendall([1.0, 2.0, 3.0], [1.0, 2.0, 3.0])
1.0

julia> kendall([1.0, 2.0, 3.0], [3.0, 2.0, 1.0])
-1.0
```


**See Also**

[`spearman`](/reference/api-metrics#Kora.spearman), [`pearson`](/reference/api-metrics#Kora.pearson)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/Kora.jl/blob/9ccdfd81b778a9c3de77df920a4498af80a8832a/src/metrics.jl#L214-L270" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Kora.ALL_METRICS' href='#Kora.ALL_METRICS'><span class="jlbinding">Kora.ALL_METRICS</span></a> <Badge type="info" class="jlObjectType jlConstant" text="Constant" /></summary>



```julia
ALL_METRICS
```


Ordered vector of all five performance metric functions used to evaluate fitted polynomial growth and survival models against EcoRRAP held-out field data.

Contains, in order: [`RMSE`](/reference/api-metrics#Kora.RMSE), [`R2`](/reference/api-metrics#Kora.R2), [`pearson`](/reference/api-metrics#Kora.pearson), [`spearman`](/reference/api-metrics#Kora.spearman), [`kendall`](/reference/api-metrics#Kora.kendall).

Each element is a callable with the shared signature `f(y_hat, y)`, where `y_hat` is the vector of model-predicted values and `y` is the vector of observed values. This uniform signature allows `ALL_METRICS` to be iterated directly when populating the performance table in `process_ecorrap_models`, and ensures the table columns produced by `show(::PolyGrowthModel)` and `show(::PolySurvivalModel)` are always in a consistent order.

**Examples**

```julia
julia> using Kora

julia> length(ALL_METRICS)
5

julia> all(f([1.0, 2.0, 3.0], [1.0, 2.0, 3.0]) isa Float64 for f in ALL_METRICS)
true
```


**See Also**

[`RMSE`](/reference/api-metrics#Kora.RMSE), [`R2`](/reference/api-metrics#Kora.R2), [`pearson`](/reference/api-metrics#Kora.pearson), [`spearman`](/reference/api-metrics#Kora.spearman), [`kendall`](/reference/api-metrics#Kora.kendall)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/Kora.jl/blob/9ccdfd81b778a9c3de77df920a4498af80a8832a/src/metrics.jl#L273-L303" target="_blank" rel="noreferrer">source</a></Badge>

</details>

