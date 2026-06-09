# Getting Started

This page walks through the minimum steps needed to install Kora, load a reef, run a simulation, and read the output. It assumes you have Julia 1.9 or later installed. For background on what the model does and why, see [What Can Kora Tell Me?](what-can-kora-tell-me.md) and the [Model Overview](model-overview.md).

## Installation

Install Kora from the Julia package registry using the standard package manager.

```julia
] add Kora
```

Alternatively, from within a script or the REPL:

```julia
using Pkg
Pkg.add("Kora")
```

## Loading the Package

Load Kora with `using Kora`. The bundled growth and survival models are loaded automatically at startup and are available as `Kora.growth_models` and `Kora.survival_models`. You do not need to load them manually.

These bundled models were fitted to offshore northern Great Barrier Reef sites. They are provided for demonstration purposes only. Any serious application should use models fitted to survey data from the target region. See the [Fitting Models from EcoRRAP Data](tutorials/fitting-from-ecorrap.md) tutorial for how to fit and supply your own models.

```julia
using Kora
```

## Creating a Reef

`initialize_reef` allocates a `ReefState` object that holds all simulation state. The call below sets up a reef with 50 annual timesteps, 10 locations, a reef area of 200 $\text{m}^2$ per location, and an initial colony density of 15 colonies per $\text{m}^2$.

```julia
reef = initialize_reef(;
    n_timesteps = 50,
    n_locs      = 10,
    area        = 200.0,
    density     = 15,
    growth_models   = Kora.growth_models,
    survival_models = Kora.survival_models
)
```

After `initialize_reef` returns, the population arrays are empty. Call `initialize_coral_population!` to seed the starting colonies before running the model.

```julia
initialize_coral_population!(reef)
```

## Generating Example Environmental Forcing

Kora requires DHW (degree heating week) data as environmental forcing. For testing and exploration, `generate_example_environment` generates a synthetic DHW time series with realistic warming trends and heatwave events.

```julia
environ = generate_example_environment(50, 10)
```

The first argument is the number of timesteps and the second is the number of locations. Both must match the values used in `initialize_reef`. The function returns a `YAXArray` that can be passed directly to `run_model!`.

For real DHW data, see the [Running Simulations](tutorials/running-simulations.md) tutorial.

## Running the Model

`run_model!` advances the simulation forward through all timesteps. It modifies the reef state in place. After it returns, the reef state object contains the complete simulation history.

```julia
run_model!(reef, environ)
```

## Reading the Output

Two functions cover the most common output needs.

`coral_cover(reef)` returns a `Vector{Float32}` of total coral cover in m^2, with one value per timestep summed across all locations.

```julia
cover = coral_cover(reef)
```

`group_cover(reef)` returns a `Matrix{Float32}` of shape `(n_timesteps, n_groups)` showing mean cover per functional group across all locations at every timestep.

```julia
gc = group_cover(reef)
```

For more detail on interpreting these outputs, querying cover at a specific timestep, and working with juvenile cover, see the [Running Simulations](tutorials/running-simulations.md) tutorial.

## Visualizing Results

Kora includes built-in visualization functions powered by Makie. These render publication-quality figures and can be displayed interactively or saved to disk.

To create a simple timeseries plot of total coral cover and DHW conditions:

```julia
using Kora.viz

fig = Kora.viz.timeseries(reef, environ)
display(fig)
```

To save a figure instead of displaying it:

```julia
save("coral_cover.png", fig)
```

For ensemble simulations, visualize multiple runs together with quantile bands:

```julia
# After running an ensemble
fig = Kora.viz.ensemble_timeseries(reef, results, environ)
```

Additional visualization functions are available for group-level cover, juvenile dynamics, and detailed analysis. Plotly support may be added in the future. For a complete list of visualization functions and examples, see the [API Reference](reference/api-reef-state.md).

