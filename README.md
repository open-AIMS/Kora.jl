# Kora.jl

A Julia framework for simulating coral reef population dynamics, including growth, mortality, bleaching, and larval recruitment. Designed for use with [EcoRRAP](https://apps.aims.gov.au/metadata/view/5d14f00f-6c24-43c1-b44e-70fe013d0757) field survey data to fit empirically grounded models and evaluate restoration scenarios under climate stress.

## Overview

Kora models five coral functional groups across an arbitrary reef area:

- Tabular *Acropora*
- Corymbose *Acropora*
- Branching non-*Acropora*
- Small massives
- Large massives

Each timestep (annual) applies growth, bleaching mortality, background mortality, and larval recruitment. Simulations can be run as single trajectories or as ensembles for parameter uncertainty analysis. Coral restoration via nursery-raised deployment is supported.

Key processes:

- **Growth**: Polynomial regression on empirical size-growth data, constrained by available habitat
- **Bleaching mortality**: Size- and depth-dependent, driven by Degree Heating Weeks (DHW)
- **Recruitment**: Local self-seeding plus external larval supply; genetic selection via the Breeder's equation
- **Restoration**: Deploy nursery-raised corals of specified functional groups and volumes

## Requirements

- Julia ~1.11.7

## Setup

### 1. Clone the repository

```bash
git clone https://github.com/ConnectedSystems/Kora.jl.git
cd Kora.jl
```

### 2. Instantiate the Julia environment

```julia
julia --project=.

# In the Julia REPL:
using Pkg
Pkg.instantiate()
```

This installs all dependencies listed in [Project.toml](Project.toml).

### 3. (Optional) Install visualization extensions

For plotting support, install one or both optional backends:

```julia
Pkg.add("Makie")          # 2D/3D interactive plots
Pkg.add("UnicodePlots")   # Terminal plots
Pkg.add("Term")           # Required alongside UnicodePlots
```

## Usage

### Programmatic usage

```julia
using Kora

# Initialize reef state (75 timesteps, 1 location, 60 m² area)
reef = initialize_reef(n_timesteps=75, n_locs=1, area=60.0, density=10, depths=7.0)

# Generate example environmental conditions (DHW time series)
env = generate_example_environment(75, 1)

# Seed initial coral population (equal proportions across 5 functional groups)
initialize_coral_population!(reef, 1, ceil(Int64, 3 * 60.0); group_proportions=fill(0.2f0, 5))

# Run a single simulation
run_example!(reef, env)

# Extract coral cover results
cover = coral_cover(reef)
```

### Prototype Interactive dashboard

The primary interface is a web dashboard showing coral cover trajectories over a 75-year horizon with sliders for restoration deployment parameters:

```bash
julia --project=. bin/main.jl
```

Open `http://localhost:9384` in a browser. The dashboard requires Makie (specifically `WGLMakie`) and `Bonito`.

### Fitting growth/survival models from EcoRRAP survey data

```julia
results = process_growth_models(
    "path/to/ecorrap_adult_juv_combined.csv",
    "path/to/ecorrap_to_species.csv";
    region="offshore_north",
    growth_degree=1
)

growth_models   = results.growth_fits
survival_models = results.survival_fits
```

Pre-fitted model objects can be serialized to [assets/models/](assets/models/) for reuse.

### Ensemble simulations

```julia
# 100 parameter sets across 22 parameters
ensemble_params = rand(22, 100)
run_ensemble!(reef, env, ensemble_params)
```

## Project Structure

```
Kora.jl/
├── src/
│   ├── Kora.jl          # Module entry point
│   ├── stats.jl              # Statistical utilities
│   ├── metrics.jl            # RMSE, R², correlation metrics
│   ├── corals/               # Coral biology models
│   │   ├── growth_model.jl   # Polynomial growth with habitat constraint
│   │   ├── mortality_model.jl# DHW/size/depth-dependent bleaching
│   │   ├── recruitment.jl    # Larval production and genetic selection
│   │   └── size_classes.jl   # Size distribution management
│   ├── reefs/                # Reef-level simulation
│   │   ├── ReefState.jl      # Core state data structure
│   │   ├── reef_dynamics.jl  # Annual timestep dynamics
│   │   ├── run_example.jl    # Single simulation runner
│   │   └── run_ensemble.jl   # Ensemble runner
│   └── interface/            # Data ingestion and model fitting
│       ├── observations.jl   # EcoRRAP data processing
│       ├── regressions.jl    # Growth and survival model fitting
│       └── create_models.jl  # High-level model creation API
├── ext/
│   ├── MakieExt/             # Makie visualization extension
│   └── UnicodePlotsExt/      # Terminal plotting extension
├── bin/
│   └── main.jl               # Interactive web dashboard
└── assets/
    ├── target_groups.csv     # Functional group definitions
    └── models/               # Serialized fitted models
```
