# Kora.jl

Kora is a simulation model for decision support in coral reef restoration planning.
It is designed for use under deep uncertainty: it maps how outcomes vary across climate
and ecological scenarios rather than forecasting a single future.

Crown-of-thorns, cyclones, flood plumes, disease, macroalgae competition, and local
stressors are not currently represented. Growth and survival models are extensible per
functional group; disturbances are currently implemented case-by-case and represent a
known extensibility target.

See [What Can Kora Tell Me?](what-can-kora-tell-me.md) for an overview of model scope and purpose.

## Kora and ADRIA in a Decision Workflow

Kora and ADRIA serve different roles in a shared decision-support workflow.

Kora provides the coral ecology representation. It simulates coral population dynamics and resulting cover trajectories under alternative climate, thermal adaptation, and restoration assumptions. These outputs are scenario-conditioned ecological outcomes, not forecasts of the most likely reef future.

ADRIA.jl provides the surrounding assessment layer for scenario ensembles, sensitivity analysis, and comparison of strategy performance across plausible futures. In this context, robustness means identifying options that remain acceptable, or fail less severely, across a defined set of plausible futures.

Current integration is workflow-level: Kora outputs can be transferred into ADRIA analyses, but there is not yet a unified package-level interface. An open Kora.jl issue proposes adopting data structures from open-AIMS/ADRIAIndicators.jl to reduce translation effort and improve interoperability. That alignment is planned work, not a completed interface.

- Use Kora to define ecological scenarios and simulate coral outcome trajectories.
- Use ADRIA to structure scenario ensembles, assess sensitivity to assumptions, and compare decision performance across plausible futures.
- Feed insights back into Kora scenario design for iterative refinement of ecological assumptions and intervention options.

See the [ADRIA.jl repository](https://github.com/open-AIMS/ADRIA.jl) and [ADRIAIndicators.jl](https://github.com/open-AIMS/ADRIAIndicators.jl).

## Quick Start

For Julia newcomers, the steps below provide a minimal local install and first simulation workflow.

### 1. Install Julia

Install Julia 1.11+ with `juliaup` (recommended Julia installer and version manager):

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

Paste and run:

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

To create a quick plot of the run:

```julia
using CairoMakie

fig = Kora.viz.timeseries(reef, environ)
display(fig)
```

For a more complete walk-through, see [Running Simulations](tutorials/running-simulations.md) and [Visualizing Simulation Results](tutorials/visualization.md).
