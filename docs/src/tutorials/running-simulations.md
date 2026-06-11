# Running Simulations

This tutorial covers the full simulation workflow in detail. It builds on [Getting Started](../getting-started.md), which shows the minimum steps. This page explains each argument, describes the output functions, and shows how to use real DHW data.

## Setting Up a Reef

`initialize_reef` allocates a `ReefState` and configures the simulation dimensions. All keyword arguments have defaults; supply only the arguments that need to be overridden.

```julia
using Kora

reef = initialize_reef(;
    n_timesteps=75,
    n_locs=30,
    area=250.0,
    density=18,
    depths=9.0,
    growth_models=Kora.growth_models,
    survival_models=Kora.survival_models
)
```

The arguments are described below.

`n_timesteps` sets the number of annual timesteps. The default is 75.

`n_locs` sets the number of reef locations. All locations run independently within a single simulation.

`area` is the reef area available for coral cover in m^2. Provide a scalar to apply the same value to every location, or a `Vector` with one entry per location.

`density` is the maximum colony density in colonies per m^2. This sets the carrying capacity ceiling together with `area`. Provide a scalar or a `Vector{Int64}` with one entry per location.

`depths` is the water depth in meters. Depth controls the bleaching coefficient: deeper locations experience reduced bleaching for a given DHW value. Provide a scalar to use the same depth for every location, or a `Vector{Float64}`.

`growth_models` and `survival_models` are the fitted polynomial model collections. The package-level defaults cover offshore northern Great Barrier Reef sites. Pass a collection loaded with `load_models` to use region-specific parameters.

Here is an example with per-location areas and depths.

```julia
areas  = fill(200.0, 30)
depths = rand(5.0:1.0:15.0, 30)

reef = initialize_reef(;
    n_timesteps=75,
    n_locs=30,
    area=areas,
    density=18,
    depths=depths,
    growth_models=Kora.growth_models,
    survival_models=Kora.survival_models
)
```

## Initialising the Coral Population

`initialize_reef` returns an empty reef. The starting population must be seeded before the model can run.

```julia
initialize_coral_population!(reef)
```

This draws initial colony diameters from per-group log-normal size distributions and writes them into `reef.wild_population` at timestep 1. The total colony count is derived from the reef area and density. Group proportions default to `[0.10, 0.20, 0.25, 0.20, 0.25]` across the five functional groups.

For reproducible initial populations, pass a seeded random number generator.

```julia
using Random
rng = Random.default_rng()
initialize_coral_population!(reef; rng=rng)
```

## Preparing Environment Data

### Synthetic data for testing

`generate_example_environment` generates a synthetic DHW time series with a warming trend, seasonal cycles, and acute heatwave events. It is useful for testing and exploring model behaviour before working with real data.

```julia
environ = generate_example_environment(75, 30)
```

The first argument is the number of timesteps and the second is the number of locations. Both must match the values used in `initialize_reef`. The function returns a `YAXArray` with axes `(Dim{:timestep}, Dim{:location}, Dim{:variable})`.

### Real DHW data

Real DHW data requires a `Matrix{Float32}` of shape `(n_timesteps, n_locs)` where rows are years and columns are reef locations. Then pass it to `generate_environment`.

```julia
dhw_matrix = Matrix{Float32}(site_dhw_data)   # shape: (75, 30)
environ = generate_environment(dhw_matrix)
```

`generate_environment` validates the data before wrapping it and will issue advisory warnings in two situations. If the maximum DHW value exceeds 40 deg-weeks, the warning indicates a likely unit mismatch, as that value is roughly twice the DHW projected under the most severe emissions scenarios. If the minimum DHW value across all locations and timesteps exceeds 20, the warning indicates that raw sea-surface temperature may have been passed instead of accumulated DHW values. Real DHW data always includes near-zero values during non-bleaching years and in early simulation periods.

## Running the Model

`run_model!` advances the simulation through all timesteps and writes results into the reef state object in place.

```julia
run_model!(reef, environ)
```

The `rng` keyword enables reproducible stochastic draws for survival and recruitment.

```julia
rng = Random.default_rng()
run_model!(reef, environ; rng=rng)
```

Two optional keyword arguments control recruitment behaviour. `recruits` is the fraction of local larval production that successfully recruits each timestep (default: `0.06`). `self_seed` is the fraction of recruitment attributed to self-seeding from the local population (default: `0.3`).

```julia
run_model!(reef, environ; recruits=0.05f0, self_seed=0.25f0)
```

After `run_model!` returns, the reef state contains the full time series for all locations, groups, and timesteps.

## Reading Results

The following functions extract summary output from a completed reef state.

`coral_cover(reef)` returns a `Vector{Float32}` with one value per timestep. Each value is the total coral cover in m^2, summed across all locations at that timestep.

```julia
cover = coral_cover(reef)   # length == n_timesteps
```

To get per-location cover at a specific timestep, pass the timestep index.

```julia
cover_ts10 = coral_cover(reef, 10)   # Vector{Float32} of length n_locs
```

`group_cover(reef)` returns a `Matrix{Float32}` of shape `(n_timesteps, n_groups)` with mean cover per functional group averaged across all locations.

```julia
gc = group_cover(reef)   # (n_timesteps, 5)
```

`group_cover_timeseries(reef)` returns the same matrix. It is the function that `group_cover(reef)` delegates to when no timestep is supplied.

`juvenile_cover(reef, ts)` returns a `Vector{Float32}` with one element per functional group, showing mean sub-mature coral cover in m^2 at timestep `ts`, averaged across all locations. Juveniles are colonies whose diameter is below the maturity threshold for their group.

```julia
jc = juvenile_cover(reef, 10)   # Vector of length n_groups
```

## Running Ensembles

In practice, many simulations are often needed to explore how outcomes vary across uncertain parameter combinations or environmental scenarios. The `run_ensemble!` function automates this: it accepts a parameter matrix where each column is one scenario and runs one simulation per column.

The parameter matrix has shape `(n_params, n_members)`. Rows 1–16 control the population state (density, group proportions, initial size distributions). Rows 17–21 are per-group location scalers (multipliers applied to growth at each location). Rows 22–23 are `recruits` and `self_seed` parameters.

Here is a simple example that runs 10 scenarios with different initial densities while holding everything else constant.

```julia
using Random

rng = Random.default_rng()

n_scenarios = 10

# Build parameter matrix: 23 rows (23 parameters), 10 columns (10 scenarios)
params = zeros(Float64, 23, n_scenarios)

# Vary density across scenarios (parameter 1)
params[1, :] = range(10.0, 25.0; length=n_scenarios)

# Run ensemble
results = run_ensemble!(reef, environ, params; rng=rng)

# Results is a NamedTuple containing:
#   results.cover        - shape (n_timesteps, n_locs, n_scenarios)
#   results.group_cover  - shape (n_timesteps, n_groups, n_scenarios)
#   results.juvenile_cover  - shape (n_timesteps, n_groups, n_scenarios)
```

For reproducible ensemble members, seed the RNG. Each column of the parameter matrix corresponds to one ensemble member and will be run in parallel if Julia is started with multiple threads.
