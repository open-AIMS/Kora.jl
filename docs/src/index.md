# Kora.jl

Kora is a simulation model for decision support in coral reef restoration planning.
It is designed for use under deep uncertainty: it maps how outcomes vary across climate
and ecological scenarios rather than forecasting a single future.

Crown-of-thorns, cyclones, flood plumes, disease, macroalgae competition, and local
stressors are not currently represented. Growth and survival models are extensible per
functional group; disturbances are currently implemented case-by-case and represent a
known extensibility target.

See [What Can Kora Tell Me?](what-can-kora-tell-me.md) for an overview of model scope and purpose.

## Quick Start

For Julia newcomers, the steps below provide a minimal local install and first simulation workflow.

### 1. Install Julia

Install Julia 1.12+ with `juliaup` (recommended Julia installer and version manager):

- https://julialang.org/install/

If `juliaup` is already installed, install or update Julia from a terminal:

```bash
juliaup add release
juliaup default release
```

After installation, open a terminal and start Julia:

```bash
julia
```

### 2. Install Kora

In the Julia REPL, press `]` to enter package mode. The prompt changes from `julia>` to `pkg>`.

Standard install (from Julia General registry, once Kora is registered):

```julia
add Kora
```

Current install (development version from GitHub):

```julia
add https://github.com/open-AIMS/Kora.jl
```

Press Backspace or `Ctrl+C` to leave package mode and return to the `julia>` prompt.

### 3. Run a first model

Example:

```julia
using Kora
using Random

rng = Random.default_rng()

reef = initialize_reef(; n_timesteps=75, n_locs=1, area=60.0, density=10, depths=7.0)
initialize_coral_population!(reef; rng)
environ = generate_example_environment(75, 1; rng)

run_model!(reef, environ; rng)
cover = coral_cover(reef)

println("Final coral cover = $(cover[end])")
```

### 4. Plot results

A quick plot of the run:

```julia
using CairoMakie

fig = Kora.viz.timeseries(reef, environ)
display(fig)
```

For a more complete walk-through, see [Running Simulations](tutorials/running-simulations.md) and [Visualizing Simulation Results](tutorials/visualization.md).

## Kora in a Broader Workflow

Kora is often used alongside [ADRIA.jl](https://github.com/open-AIMS/ADRIA.jl) as part of a larger decision-support pipeline. See [Kora and ADRIA in a Decision Workflow](concepts/kora-and-adria.md) for how the two packages relate.
