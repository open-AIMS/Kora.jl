
# Simulation API {#Simulation-API}
<details class='jldocstring custom-block' open>
<summary><a id='Kora.run_model' href='#Kora.run_model'><span class="jlbinding">Kora.run_model</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
run_model(;
    n_ts=75,
    n_locs=100,
    with_dhw=true,
    area=100.0,
    pop_density=15.0,
    growth_models=growth_models,
    survival_models=survival_models
)::Tuple{ReefState, YAXArray}
```


Convenience wrapper that allocates a reef, seeds its population, generates synthetic environmental conditions, and returns the completed simulation results.

Internally calls `initialize_reef`, `initialize_coral_population!`, `generate_example_environment`, and `run_model!` in sequence using the supplied keyword arguments.

**Arguments**
- `n_ts` : Number of annual time steps (default: `75`).
  
- `n_locs` : Number of reef locations (default: `100`).
  
- `with_dhw` : Whether to generate DHW thermal forcing. Pass `false` to run with zero thermal stress (default: `true`).
  
- `area` : Reef area in m^2 used for carrying capacity (default: `100.0`).
  
- `pop_density` : Initial colony density in colonies per m^2 used to size the starting population (default: `15.0`).
  
- `growth_models` : Fitted growth model collection (default: package-level offshore-north models).
  
- `survival_models` : Fitted survival model collection (default: package-level offshore-north models).
  

**Returns**

`Tuple{ReefState, YAXArray}` : The completed reef state containing the full time series and the environmental conditions used for the run.

**See Also**

[`run_model!`](/reference/api-simulation#Kora.run_model!), [`initialize_reef`](/reference/api-reef-state#Kora.initialize_reef), [`generate_example_environment`](/reference/api-reef-state#Kora.generate_example_environment)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/Kora.jl/blob/319707aa5d4659465eb42066e8cff7c4b1e14e2e/src/reefs/run_model.jl#L297-L335" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Kora.run_model!' href='#Kora.run_model!'><span class="jlbinding">Kora.run_model!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
run_model!(
    reef_state::ReefState,
    env_conditions::YAXArray;
    recruits=0.06f0,
    self_seed=0.3f0,
    rng::AbstractRNG=Random.GLOBAL_RNG
)::Nothing
```


Advance the reef simulation forward in time, writing results into `reef_state` in-place. The state is reset to its timestep-1 population before the run begins, so calling `run_model!` a second time on the same object produces a fresh result.

`env_conditions` must be a 3D `YAXArray` with axes `(Dim{:timestep}, Dim{:location}, Dim{:variable})` containing at minimum a `:dhw` variable slice. This is the format returned by both `generate_example_environment` and `generate_environment`.

**Arguments**
- `reef_state` : Pre-initialised `ReefState`. Modified in-place; the caller's object contains the full time series after this call returns.
  
- `env_conditions` : Environmental forcing data. Must cover the same number of timesteps and locations as `reef_state`.
  
- `recruits` : Fraction of local larval production that successfully recruits to the reef each timestep (default: `0.06`).
  
- `self_seed` : Fraction of recruitment attributed to self-seeding from the local population (default: `0.3`).
  
- `rng` : Random number generator. Pass a seeded `Xoshiro` or similar for reproducible runs (default: `Random.GLOBAL_RNG`).
  

**Returns**

`Nothing`

**See Also**

[`initialize_reef`](/reference/api-reef-state#Kora.initialize_reef), [`initialize_coral_population!`](/reference/api-reef-state#Kora.initialize_coral_population!), [`run_model`](/reference/api-simulation#Kora.run_model), [`coral_cover`](/reference/api-reef-state#Kora.coral_cover)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/Kora.jl/blob/319707aa5d4659465eb42066e8cff7c4b1e14e2e/src/reefs/run_model.jl#L1-L37" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Kora.run_ensemble!' href='#Kora.run_ensemble!'><span class="jlbinding">Kora.run_ensemble!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
run_ensemble!(
    reef_state::ReefState,
    env_conditions::YAXArray,
    ensemble_params::Matrix{Float64};
    rng::AbstractRNG=Random.GLOBAL_RNG
)
```


Run an ensemble of simulations, one per column of `ensemble_params`, reusing `reef_state` across members (reset between each run via `set_population!`). Progress is reported to the terminal via `ProgressLogging`.

When `ensemble_params` has more than 16 rows, rows 17 through `16 + n_groups(reef_state)` are interpreted as per-group location scalers passed to `assign_scalers!`, and the final two rows as `recruits` and `self_seed` parameters forwarded to `run_model!`.

**Arguments**
- `reef_state` : `ReefState` used as the simulation template. Mutated during each member run; contents after the call reflect only the last ensemble member.
  
- `env_conditions` : Environmental forcing data shared across all members.
  
- `ensemble_params` : Parameter matrix of shape `(n_params, n_members)`. Each column defines one ensemble member. Rows 1-16 are population parameters consumed by `set_population!`.
  
- `rng` : Random number generator (default: `Random.GLOBAL_RNG`).
  

**Returns**

`NamedTuple` with the fields listed below.
- `cover::Array{Float32,3}` : Total coral cover in m^2 with shape `(n_timesteps, n_locations, n_members)`.
  
- `group_cover::Array{Float32,4}` : Per-group cover in m^2 with shape `(n_timesteps, n_locations, n_groups, n_members)`.
  
- `juvenile_cover::Array{Float32,4}` : Sub-mature coral cover in m^2 with shape `(n_timesteps, n_locations, n_groups, n_members)`.
  
- `wild_dhw_tolerances::Array{Float32,5}` : DHW tolerance statistics with shape `(n_timesteps, n_locations, n_groups, 2, n_members)`. The third inner dimension holds mean (index 1) and standard deviation (index 2).
  
- `params::Matrix{Float64}` : The input `ensemble_params` unchanged.
  

**See Also**

[`run_model!`](/reference/api-simulation#Kora.run_model!), [`set_population!`](/reference/api-simulation#Kora.set_population!), [`assign_scalers!`](/reference/api-simulation#Kora.assign_scalers!)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/Kora.jl/blob/319707aa5d4659465eb42066e8cff7c4b1e14e2e/src/reefs/run_ensemble.jl#L3-L45" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Kora.set_population!' href='#Kora.set_population!'><span class="jlbinding">Kora.set_population!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
set_population!(reef_state::ReefState, x::Vector)::Nothing
```


Set the initial population state.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/Kora.jl/blob/319707aa5d4659465eb42066e8cff7c4b1e14e2e/src/reefs/interface.jl#L24-L28" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Kora.assign_scalers!' href='#Kora.assign_scalers!'><span class="jlbinding">Kora.assign_scalers!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
assign_scalers!(reef_state::ReefState, x::Vector)::Nothing
```


Assign growth and survival scalers, assuming they are identical for each location.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/Kora.jl/blob/319707aa5d4659465eb42066e8cff7c4b1e14e2e/src/reefs/interface.jl#L3-L7" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Deployment (Internal) {#Deployment-Internal}

`deploy_corals!` is called internally by `run_model!` when coral outplanting is active. It is documented here for completeness; users should configure deployment through the `run_model!` interface rather than calling it directly.
<details class='jldocstring custom-block' open>
<summary><a id='Kora.deploy_corals!' href='#Kora.deploy_corals!'><span class="jlbinding">Kora.deploy_corals!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
deploy_corals!(reef_state, ts, loc, n, grp; rng=Random.GLOBAL_RNG)
```


Seed `n` outplanted coral colonies of functional group `grp` at location `loc` and timestep `ts`.

Colony initial diameters are sampled from the truncated log-normal size distribution for `grp` (bounded by the group's diameter bin edges) and stored in `reef_state.deployed_population[ts, loc, grp]`. Any existing deployed population at that slot is overwritten.

This function is called internally when coral deployment is active in `run_model!`. Direct calls are supported for testing but are not part of the standard simulation workflow – use the deployment configuration in `run_model!` to trigger outplanting within a simulation run.

**Arguments**
- `reef_state` : `ReefState` to update in-place.
  
- `ts` : Timestep index at which deployment occurs (1-based).
  
- `loc` : Location index at which corals are deployed (1-based).
  
- `n` : Number of colonies to deploy.
  
- `grp` : Functional group index (1-based).
  
- `rng` : Random number generator for reproducible diameter draws (default: `Random.GLOBAL_RNG`).
  

**Returns**

`Nothing`

**See Also**

[`initialize_coral_population!`](/reference/api-reef-state#Kora.initialize_coral_population!), [`run_model!`](/reference/api-simulation#Kora.run_model!)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/Kora.jl/blob/319707aa5d4659465eb42066e8cff7c4b1e14e2e/src/reefs/ReefState.jl#L567-L597" target="_blank" rel="noreferrer">source</a></Badge>

</details>

