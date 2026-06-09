
# Reef State API {#Reef-State-API}
<details class='jldocstring custom-block' open>
<summary><a id='Kora.ReefState' href='#Kora.ReefState'><span class="jlbinding">Kora.ReefState</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
ReefState{F, P, Y3a, Y4a, Y4b}
```


Mutable container holding the full ecological state of a simulated reef system across time, space, and functional coral groups.

**Fields**
- `wild_population::Array{Vector{F},3}` : Diameter samples (cm) for wild corals, indexed `[timestep, location, group]`. Each element is a `Vector` of individual colony diameters at that point in the simulation.
  
- `deployed_population::Array{Vector{F},3}` : Diameter samples (cm) for outplanted corals, with the same `[timestep, location, group]` indexing.
  
- `deployment_times::Array{F,3}` : Number of corals deployed at each `[timestep, location, group]` combination.
  
- `growth_models::Vector{P}` : Per-group growth functions. Each callable maps a colony diameter (cm) to its expected diameter at the next annual timestep.
  
- `survival_models::Vector{P}` : Per-group survival functions. Each callable maps a colony diameter (cm) to an annual survival probability.
  
- `carrying_capacity::Vector{F}` : Maximum coral-bearing area in m^2 for each location. Limits total cover and constrains recruitment.
  
- `depths::Vector{F}` : Water depth in meters at each location. Used to compute depth-dependent bleaching mortality coefficients.
  
- `density::Vector{Int64}` : Maximum colony density in colonies per m^2 for each location. Recruitment is suppressed when the total population approaches this ceiling.
  

All other fields (YAXArray tolerance and mortality stores, and fields prefixed with `_`) are internal implementation details that may change between minor versions.

**See Also**

[`initialize_reef`](/reference/api-reef-state#Kora.initialize_reef), [`run_model!`](/reference/api-simulation#Kora.run_model!)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/Kora.jl/blob/319707aa5d4659465eb42066e8cff7c4b1e14e2e/src/reefs/ReefState.jl#L5-L37" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Kora.initialize_reef' href='#Kora.initialize_reef'><span class="jlbinding">Kora.initialize_reef</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
initialize_reef(;
    n_timesteps::Int=75,
    n_locs::Int=100,
    group_names=Kora.TARGET_GROUPS,
    density::Union{Int64,Vector{Int64}}=20,
    area=90.0,
    depths::Union{Float64,Vector{Float64}}=9.0,
    growth_models::AbstractCoralBehavior=Kora.growth_models,
    survival_models::AbstractCoralBehavior=Kora.survival_models
)::ReefState
```


Allocate and return a `ReefState` sized for `n_locs` reef locations and `n_timesteps` annual time steps.

All population arrays are initialised empty. Call `initialize_coral_population!` to seed the starting population before running a simulation.

**Arguments**
- `n_timesteps` : Number of annual time steps to allocate (default: `75`).
  
- `n_locs` : Number of reef locations (default: `100`).
  
- `group_names` : Labels for the functional coral groups. Must match the groups used to fit `growth_models` and `survival_models` (default: `Kora.TARGET_GROUPS`, the five groups used by the bundled models).
  
- `density` : Maximum colony density in colonies per m^2. Provide a scalar to apply the same ceiling to every location, or a per-location `Vector{Int64}` (default: `20`).
  
- `area` : Reef area available for coral cover in m^2. Provide a scalar or a per-location `Vector` (default: `90.0`).
  
- `depths` : Water depth in meters. Provide a scalar or a per-location `Vector{Float64}`. Depth controls bleaching mortality coefficients (default: `9.0`).
  
- `growth_models` : Fitted growth model collection, one function per functional group. Defaults to the package-level offshore-north models loaded from the bundled JSON asset at package load time.
  
- `survival_models` : Fitted survival model collection, one function per functional group. Defaults to the package-level offshore-north models loaded from the bundled JSON asset at package load time.
  

**Returns**

`ReefState` : An empty reef state ready for population initialisation.

**Examples**

```julia
julia> using Kora

julia> rs = initialize_reef(; n_timesteps=10, n_locs=3);

julia> n_timesteps(rs), n_locations(rs)
(10, 3)
```


**See Also**

[`initialize_coral_population!`](/reference/api-reef-state#Kora.initialize_coral_population!), [`run_model!`](/reference/api-simulation#Kora.run_model!), [`load_models`](/reference/api-interface#Kora.load_models)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/Kora.jl/blob/319707aa5d4659465eb42066e8cff7c4b1e14e2e/src/reefs/ReefState.jl#L299-L354" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Kora.initialize_coral_population!' href='#Kora.initialize_coral_population!'><span class="jlbinding">Kora.initialize_coral_population!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
initialize_coral_population!(
    reef_state::ReefState,
    loc::Int64,
    target_pop_size::Int64;
    group_proportions::Vector{Float32}=[0.1f0, 0.2f0, 0.25f0, 0.2f0, 0.25f0],
    rng::AbstractRNG=Random.GLOBAL_RNG
)::Nothing

initialize_coral_population!(
    reef_state::ReefState;
    rng::AbstractRNG=Random.GLOBAL_RNG
)
```


Seed the initial coral population at timestep 1 with colony diameter samples drawn from per-group log-normal size distributions.

The no-location convenience method seeds all locations using a target population size derived from the maximum carrying capacity in `reef_state`. Colony sizes are drawn from group-specific truncated log-normal distributions and written as diameter vectors (cm) into `wild_population[1, loc, grp]`.

**Arguments**
- `reef_state` : The `ReefState` to populate. Population data are written in-place at timestep 1.
  
- `loc` : Index of the location to seed (1-based).
  
- `target_pop_size` : Total number of colonies to place at this location. Actual per-group counts are `round(target_pop_size * group_proportions[grp])`.
  
- `group_proportions` : Proportion of the total population assigned to each functional group. Must sum to 1.0 (checked with `atol=1e-6`). Default: `[0.10, 0.20, 0.25, 0.20, 0.25]` matching the five bundled groups.
  
- `rng` : Random number generator for reproducible diameter draws (default: `Random.GLOBAL_RNG`).
  

**Returns**

`Nothing`

**See Also**

[`initialize_reef`](/reference/api-reef-state#Kora.initialize_reef), [`run_model!`](/reference/api-simulation#Kora.run_model!)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/Kora.jl/blob/319707aa5d4659465eb42066e8cff7c4b1e14e2e/src/reefs/ReefState.jl#L478-L517" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Kora.generate_example_environment' href='#Kora.generate_example_environment'><span class="jlbinding">Kora.generate_example_environment</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
generate_example_environment(
    n_years::Int64,
    n_locations::Int64;
    rng::AbstractRNG=Random.GLOBAL_RNG,
    start_year::Int64=2020,
    with_dhw=true,
    warming_rate::Float32=0.15f0,
    seasonal_amplitude::Float32=1.2f0,
    dhw_threshold::Float32=4.0f0,
    noise_amplitude::Float32=0.9f0
)
```


Generate synthetic environmental data for coral reef modeling, specifically Degree Heating Weeks (DHW) trajectories that simulate realistic marine heatwave patterns under climate change scenarios.

**Extended help**

Generate plausible DHW time series by combining multiple environmental components:
1. **Long-term warming trend**: Simulates gradual ocean warming (like RCP4.5 scenarios)
  
2. **Seasonal cycles**: Models natural temperature variations throughout the year
  
3. **Weather noise**: Adds realistic short-term temperature fluctuations
  
4. **Spatial variation**: Creates location-specific temperature offsets
  
5. **Acute heatwave events**: Superimposes extreme marine heatwave events
  
6. **DHW accumulation**: Converts temperature anomalies to ecologically-relevant DHW values
  

The DHW calculation follows coral bleaching research where:
- DHW accumulates when temperatures exceed a threshold (default 4°C above baseline)
  
- Values decay over time when temperatures drop
  
- Extreme events can cause rapid DHW spikes that lead to mass bleaching
  

**Temperature Anomaly Construction**

For each location and timestep, the temperature anomaly is built from:

```julia
temp_anomaly = warming_trend + seasonal_cycle + weather_noise + spatial_offset
```


Where:
- **warming_trend**: Linear increase over time (`warming_rate * years_elapsed`)
  
- **seasonal_cycle**: Sinusoidal pattern with amplitude that increases over time
  
- **weather_noise**: Random normal variations that get larger in later years
  
- **spatial_offset**: Location-specific random offset (0 - 0.8°C)
  

**DHW Accumulation Rules**
- **Above threshold**: DHW accumulates as `(temp_anomaly - threshold) / 4.0`
  
- **Below threshold**: DHW decays rapidly (`previous_DHW * 0.7`)
  
- **During heating**: DHW decays slowly (`previous_DHW * 0.92`)
  
- **Soft cap**: Values above 20 DHW are dampened but can still fluctuate
  

**Acute Event Generation**

The function adds realistic extreme heatwave events:
- **Frequency**: ~1 event per year on average (`n_timesteps / 12`)
  
- **Probability**: Increases over time (simulating worsening climate)
  
- **Duration**: 2-5+ weeks, longer in later years
  
- **Intensity**: 8-25+ DHW, stronger in later years
  
- **Shape**: Rapid onset (30% of duration) → peak → gradual decline (30% of duration)
  

**Ecological Realism**

The parameters are tuned based on coral bleaching research:
- **4 DHW**: Threshold where bleaching typically begins
  
- **8+ DHW**: Significant bleaching and mortality expected
  
- **20+ DHW**: Severe bleaching events (soft-capped with fluctuations)
  

**Arguments**
- `n_years`: Number of time steps to simulate (in years)
  
- `n_locations`: Number of spatial locations across the reef system
  
- `rng`: Random number generator for reproducible results
  
- `start_year`: Starting year for the simulation (default: 2020)
  
- `with_dhw`: Whether to generate DHW data or return zeros (default: true)
  
- `warming_rate`: Rate of long-term warming per year in °C (default: 0.15°C/year)
  
- `seasonal_amplitude`: Strength of seasonal temperature cycles (default: 1.2°C)
  
- `dhw_threshold`: Temperature threshold above which DHW accumulates (default: 4.0°C)
  
- `noise_amplitude`: Magnitude of random weather variations (default: 0.9°C)
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/Kora.jl/blob/319707aa5d4659465eb42066e8cff7c4b1e14e2e/src/reefs/ReefState.jl#L985-L1061" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Kora.generate_environment' href='#Kora.generate_environment'><span class="jlbinding">Kora.generate_environment</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
generate_environment(dhw::Matrix{Float32}; start_year::Int64=2020)::YAXArray
generate_environment(dhw::YAXArray; start_year::Int64=2020)::YAXArray
```


Wrap user-supplied Degree Heating Weeks data in a correctly structured environment YAXArray, suitable for direct use with Kora.jl model runs.

`n_years` and `n_locations` are inferred from the input dimensions. Structure creation is delegated to `generate_example_environment` so dimension names and axis labels are always consistent. The synthetic DHW values produced by that function are then replaced with the caller's real data.

**Arguments**
- `dhw` : DHW data with shape `(n_timesteps, n_locations)`. Each row is one year (or timestep) and each column is one reef location. Accepts either a `Matrix{Float32}` or a 2D `YAXArray` whose first dimension is timestep and second dimension is location.
  
- `start_year` : First year label for the timestep axis (default: 2020). Passed through to `generate_example_environment` for axis labelling only; it does not alter the data.
  

**Returns**

A 3D `YAXArray` with axes `(Dim{:timestep}, Dim{:location}, Dim{:variable})` identical in structure to the output of `generate_example_environment`, with the `:dhw` variable populated from `dhw`.

**Notes**

Two advisory warnings are issued (not errors). A minimum-floor check fires when `minimum(dhw) > 20`, because real DHW data always contains near-zero values during non-bleaching periods – a uniformly high floor is the primary indicator that raw sea-surface temperature (~25-32 degrees C) was passed instead of DHW. A ceiling check fires when `maximum(dhw) > 40`, approximately twice the ~20 DHW projected under SSP5-8.5; values above this threshold are likely a data quality issue rather than an intentional scenario.

**Examples**

```julia
julia> using Kora

julia> dhw = zeros(Float32, 10, 5);

julia> env = generate_environment(dhw);

julia> size(env)
(10, 5, 1)
```


**See Also**

[`generate_example_environment`](/reference/api-reef-state#Kora.generate_example_environment)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/Kora.jl/blob/319707aa5d4659465eb42066e8cff7c4b1e14e2e/src/reefs/ReefState.jl#L1212-L1261" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Kora.coral_cover' href='#Kora.coral_cover'><span class="jlbinding">Kora.coral_cover</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
coral_cover(reef_state::ReefState)::Vector{Float32}
coral_cover(reef_state::ReefState, ts::Int64)::Vector{Float32}
coral_cover(reef_state::ReefState, ts::Int64, loc::Int64)::Float32
coral_cover(diams::AbstractVector{<:AbstractFloat})
```


Compute total coral cover in m^2, summing over all wild and deployed colonies across all functional groups.

The no-timestep form returns a `Vector` of summed cover across all locations for each timestep. The single-timestep form returns a per-location `Vector` at `ts`. The two-argument form returns a scalar for one timestep and one location. The bare-vector form computes cover from a diameter vector (cm) directly, without requiring a `ReefState`.

Colony area is computed as pi/4 * (d/100)^2 (m^2) for each diameter d in cm.

**Arguments**
- `reef_state` : Source of population data.
  
- `ts` : Timestep index (1-based).
  
- `loc` : Location index (1-based).
  
- `diams` : Vector of colony diameters in cm.
  

**Returns**

`Vector{Float32}` or `Float32` : Coral cover in m^2.

**See Also**

[`group_cover`](/reference/api-reef-state#Kora.group_cover), [`juvenile_cover`](/reference/api-reef-state#Kora.juvenile_cover), [`cover_cm_to_m2`](/reference/api-reef-state#Kora.cover_cm_to_m2)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/Kora.jl/blob/319707aa5d4659465eb42066e8cff7c4b1e14e2e/src/reefs/ReefState.jl#L672-L700" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Kora.cover_cm_to_m2' href='#Kora.cover_cm_to_m2'><span class="jlbinding">Kora.cover_cm_to_m2</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



Obtain the sum of all area (in m^2) for a given set of diameters.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/Kora.jl/blob/319707aa5d4659465eb42066e8cff7c4b1e14e2e/src/corals/size_classes.jl#L115-L117" target="_blank" rel="noreferrer">source</a></Badge>



Convert centimeter diameter to meters area (m^2)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/Kora.jl/blob/319707aa5d4659465eb42066e8cff7c4b1e14e2e/src/corals/size_classes.jl#L129-L131" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Kora.group_cover' href='#Kora.group_cover'><span class="jlbinding">Kora.group_cover</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
group_cover(reef_state::ReefState)::Matrix{Float32}
group_cover(reef_state::ReefState, ts::Int64)::Vector{Float32}
```


Compute mean coral cover in m^2 per functional group, averaged across all locations.

The no-timestep form delegates to `group_cover_timeseries` and returns results for every timestep as a matrix. The single-timestep form returns a per-group vector at `ts`.

**Arguments**
- `reef_state` : Source of population data.
  
- `ts` : Timestep index (1-based).
  

**Returns**

`Matrix{Float32}` with shape `(n_timesteps, n_groups)`, or `Vector{Float32}` with one element per functional group. Values are mean cover in m^2 across all locations.

**See Also**

[`coral_cover`](/reference/api-reef-state#Kora.coral_cover), [`group_cover_timeseries`](/reference/api-reef-state#Kora.group_cover_timeseries), [`juvenile_cover`](/reference/api-reef-state#Kora.juvenile_cover)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/Kora.jl/blob/319707aa5d4659465eb42066e8cff7c4b1e14e2e/src/reefs/ReefState.jl#L730-L753" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Kora.group_cover_timeseries' href='#Kora.group_cover_timeseries'><span class="jlbinding">Kora.group_cover_timeseries</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
group_cover_timeseries(reef_state::ReefState)::Matrix{Float32}
```


Compute mean coral cover in m^2 per functional group for every timestep, averaged across all locations.

**Arguments**
- `reef_state` : Source of population data.
  

**Returns**

`Matrix{Float32}` : Cover matrix of shape `(n_timesteps, n_groups)`. Rows are timesteps; columns are functional groups in the order defined by `reef_state.wild_population`.

**See Also**

[`group_cover`](/reference/api-reef-state#Kora.group_cover), [`coral_cover`](/reference/api-reef-state#Kora.coral_cover)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/Kora.jl/blob/319707aa5d4659465eb42066e8cff7c4b1e14e2e/src/reefs/ReefState.jl#L774-L790" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Kora.juvenile_cover' href='#Kora.juvenile_cover'><span class="jlbinding">Kora.juvenile_cover</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
juvenile_cover(
    reef_state::ReefState,
    ts::Int64;
    juvenile_threshold::Union{Nothing,Float32,Vector{Float32}}=nothing
)::Vector{Float32}
```


Compute mean cover of sub-mature corals in m^2 per functional group at a single timestep, averaged across all locations.

A colony is classified as juvenile when its diameter is strictly less than the maturity threshold for its group. The default thresholds come from `Kora.mature_size_thresholds()`; pass a custom value to override.

**Arguments**
- `reef_state` : Source of population data.
  
- `ts` : Timestep index (1-based).
  
- `juvenile_threshold` : Diameter threshold in cm below which a colony counts as juvenile. Provide a scalar `Float32` to apply the same value to all groups, a `Vector{Float32}` for per-group thresholds, or `nothing` to use the package defaults (default: `nothing`).
  

**Returns**

`Vector{Float32}` : Mean juvenile cover in m^2, one element per functional group.

**See Also**

[`group_cover`](/reference/api-reef-state#Kora.group_cover), [`coral_cover`](/reference/api-reef-state#Kora.coral_cover)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/Kora.jl/blob/319707aa5d4659465eb42066e8cff7c4b1e14e2e/src/reefs/ReefState.jl#L896-L924" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Kora.recruit_cover' href='#Kora.recruit_cover'><span class="jlbinding">Kora.recruit_cover</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
recruit_cover(recruits::Array{Float32})::Vector{Float32}
recruit_cover(recruits::Matrix{Vector{Float32}})::Vector{Float32}
recruit_cover(ecostate::ReefState, recruits::Array{Float32})
recruit_cover(ecostate::ReefState, recruits::Matrix{Vector{Float32}})
```


Compute total coral cover in m^2 for a cohort of new recruits, summed per location.

The single-argument forms allocate a fresh output vector. The two-argument forms write into `ecostate._recruit_buffer` and return that view, avoiding allocation inside the simulation loop.

Colony area is computed as `pi/4 * (d/100)^2` (m^2) for each recruit diameter `d` in cm via [`cover_cm_to_m2`](/reference/api-reef-state#Kora.cover_cm_to_m2).

**Arguments**
- `ecostate` : `ReefState` whose `_recruit_buffer` is used to store results.
  
- `recruits` : Per-location recruit diameter data. Either a 3-D `Array{Float32}` with axes `[location, group, colony]`, or a `Matrix{Vector{Float32}}` with dimensions `[location, group]`.
  

**Returns**

`Vector{Float32}` : Total recruit cover in m^2, one element per location.

**See Also**

[`coral_cover`](/reference/api-reef-state#Kora.coral_cover), [`cover_cm_to_m2`](/reference/api-reef-state#Kora.cover_cm_to_m2)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/Kora.jl/blob/319707aa5d4659465eb42066e8cff7c4b1e14e2e/src/reefs/ReefState.jl#L828-L854" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Kora.area_to_diam' href='#Kora.area_to_diam'><span class="jlbinding">Kora.area_to_diam</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
area_to_diam(area::AbstractFloat)::AbstractFloat
area_to_diam(area::AbstractString)::Union{AbstractFloat, Missing}
area_to_diam(_::Missing)::Missing
```


Convert a coral planar area (cm^2) to an equivalent circle diameter (cm).

Assumes the coral footprint is circular. The conversion solves for `d` in `area = pi * (d/2)^2`, giving `d = sqrt(4 * area / pi)`.

String inputs are parsed to `Float64`; strings that cannot be parsed are returned as `missing`.

**Arguments**
- `area`: Planar area of the coral colony in cm^2. Accepts `AbstractFloat`, `AbstractString`, or `Missing`.
  

**Returns**
- `AbstractFloat`: Equivalent circle diameter in cm when input is numeric.
  
- `Missing`: Returned when `area` is `missing` or is a non-numeric string.
  

**Examples**

```julia
julia> using Kora

julia> area_to_diam(Float64(pi))
2.0

julia> area_to_diam(missing)
missing
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/Kora.jl/blob/319707aa5d4659465eb42066e8cff7c4b1e14e2e/src/interface/observations.jl#L11-L42" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Statistical Utilities {#Statistical-Utilities}
<details class='jldocstring custom-block' open>
<summary><a id='Kora.truncated_standard_normal_mean' href='#Kora.truncated_standard_normal_mean'><span class="jlbinding">Kora.truncated_standard_normal_mean</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
truncated_standard_normal_mean(lb::F, ub::F)::F where {F<:Real}
```


Compute the mean of the standard normal distribution truncated to the interval [`lb`, `ub`].

Implementation follows Distributions.jl, excluding unused error checks. When `lb` &gt; `ub`, `ub` is returned to avoid NaN propagation.

**Arguments**
- `lb` : Lower bound of the truncated distribution.
  
- `ub` : Upper bound of the truncated distribution.
  

**Returns**

`F` : Mean of the truncated standard normal distribution, or `ub` if `lb` &gt; `ub`.

**Examples**

```julia
julia> using Kora

julia> truncated_standard_normal_mean(-1.0, 1.0)
0.0

julia> truncated_standard_normal_mean(0.0, 0.0)
0.0
```


**References**
1. Distributions.jl truncated normal implementation: https://github.com/JuliaStats/Distributions.jl/blob/c1705a3015d438f7e841e82ef5148224813831e8/src/truncated/normal.jl#L24-L46
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/Kora.jl/blob/319707aa5d4659465eb42066e8cff7c4b1e14e2e/src/stats.jl#L61-L91" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Kora.truncated_normal_mean' href='#Kora.truncated_normal_mean'><span class="jlbinding">Kora.truncated_normal_mean</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
truncated_normal_mean(
    normal_mean::F,
    normal_stdev::F,
    lower_bound::F,
    upper_bound::F
)::F where {F<:AbstractFloat}
```


Compute the mean of the normal distribution with mean `normal_mean` and standard deviation `normal_stdev`, truncated to the interval [`lower_bound`, `upper_bound`].

Delegates to [`truncated_standard_normal_mean`](/reference/api-reef-state#Kora.truncated_standard_normal_mean) after standardising the bounds to the unit-normal scale.

**Arguments**
- `normal_mean` : Mean of the underlying (untruncated) normal distribution.
  
- `normal_stdev` : Standard deviation of the underlying (untruncated) normal distribution.
  
- `lower_bound` : Lower bound of the truncated normal distribution.
  
- `upper_bound` : Upper bound of the truncated normal distribution.
  

**Returns**

`F` : Mean of the truncated normal distribution.

**Examples**

```julia
julia> using Kora

julia> truncated_normal_mean(0.0, 1.0, -1.0, 1.0)
0.0

julia> truncated_normal_mean(5.0, 2.0, 5.0, 5.0)
5.0
```


**See Also**

[`truncated_standard_normal_mean`](/reference/api-reef-state#Kora.truncated_standard_normal_mean), [`truncated_normal_cdf`](/reference/api-reef-state#Kora.truncated_normal_cdf)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/Kora.jl/blob/319707aa5d4659465eb42066e8cff7c4b1e14e2e/src/stats.jl#L116-L152" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Kora.truncated_normal_cdf' href='#Kora.truncated_normal_cdf'><span class="jlbinding">Kora.truncated_normal_cdf</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
truncated_normal_cdf(
    x::F,
    normal_mean::F,
    normal_stdev::F,
    lower_bound::F,
    upper_bound::F
)::F where {F<:Real}
```


Evaluate the CDF of the normal distribution with mean `normal_mean` and standard deviation `normal_stdev`, truncated to the interval [`lower_bound`, `upper_bound`], at the point `x`.

Returns `0.0` when `x <= lower_bound` and `1.0` when `x >= upper_bound`. Uses a rational approximation of the error function via `rational_erf` for efficiency; falls back to `SpecialFunctions.erf` when the truncation bounds lie more than 3 standard deviations from the mean to avoid precision loss.

**Arguments**
- `x` : Value at which to evaluate the CDF.
  
- `normal_mean` : Mean of the underlying (untruncated) normal distribution.
  
- `normal_stdev` : Standard deviation of the underlying (untruncated) normal distribution.
  
- `lower_bound` : Lower bound of the truncated distribution.
  
- `upper_bound` : Upper bound of the truncated distribution.
  

**Returns**

`F` : CDF of the truncated normal distribution evaluated at `x`, in [0, 1].

**Examples**

```julia
julia> using Kora

julia> truncated_normal_cdf(-2.0, 0.0, 1.0, -1.0, 1.0)
0.0

julia> truncated_normal_cdf(2.0, 0.0, 1.0, -1.0, 1.0)
1.0

julia> truncated_normal_cdf(0.0, 0.0, 1.0, -1.0, 1.0)
0.5
```


**See Also**

[`truncated_normal_mean`](/reference/api-reef-state#Kora.truncated_normal_mean), [`truncated_standard_normal_mean`](/reference/api-reef-state#Kora.truncated_standard_normal_mean)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/Kora.jl/blob/319707aa5d4659465eb42066e8cff7c4b1e14e2e/src/stats.jl#L162-L206" target="_blank" rel="noreferrer">source</a></Badge>

</details>

