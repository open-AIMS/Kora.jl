# Kora.jl

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.20300715-blue.svg)](https://doi.org/10.5281/zenodo.20300715)
[![Documentation (dev)](https://img.shields.io/badge/docs-dev-blue.svg)](https://open-aims.github.io/Kora.jl/dev/)

A Julia framework for simulating coral reef population dynamics, including growth, mortality, bleaching, and larval recruitment. Designed for decision support purposes, with use of [EcoRRAP](https://apps.aims.gov.au/metadata/view/5d14f00f-6c24-43c1-b44e-70fe013d0757) field survey data to fit empirically grounded models and evaluate restoration scenarios under climate stress.

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

- Julia ~1.11.7 or 1.12+ (1.12 required for AOT compilation)

## Development Setup

### 1. Clone the repository

```bash
git clone https://github.com/ConnectedSystems/Kora.jl.git
cd Kora.jl
```

### 2. Instantiate the Julia environment

```julia
julia --project=.

# In the Julia REPL:
] instantiate
```

This installs all dependencies listed in [Project.toml](Project.toml).

### 3. (Optional) Install visualization extensions

For plotting support, install one or both optional backends:

```julia
] add CairoMakie  # One of the Makie backends

# Or for plots in terminal
] add UnicodePlots Term
```

## AOT Compilation (kora_ui)

`bridge_aot.jl` and `build_julia.jl` support compiling Kora into a native shared library for use by [kora_ui](https://github.com/ConnectedSystems/kora_ui), the Rust/egui frontend. The output (`libkora_bridge.dll/.so/.dylib`) exposes a C ABI callable without the Julia runtime startup cost.

**Requires Julia 1.12+** (`juliac` must be on `PATH` or discoverable via `Sys.BINDIR`).

```powershell
# Windows (PowerShell) — from the Kora.jl root:
.\build\build.ps1

# Override output directory:
.\build\build.ps1 -OutputDir C:\path\to\kora_ui\julia_lib

# Or via environment variable:
$env:KORA_LIB_DIR = "C:\path\to\kora_ui\julia_lib"; .\build\build.ps1
```

```bash
# Linux / macOS — from the Kora.jl root:
./build/build.sh

# Override output directory:
./build/build.sh --output-dir /path/to/kora_ui/julia_lib

# Or via environment variable:
KORA_LIB_DIR=/path/to/kora_ui/julia_lib ./build/build.sh
```

Output is written to `julia_lib/` by default. The script invokes:

```
juliac --project=. --output-lib julia_lib/kora_bridge --trim=safe --compile-ccallable --experimental build/bridge_aot.jl
```

### C API

`bridge_aot.jl` exports two `@ccallable` functions. State is global; `kf_init_reef` must be called before `kf_run_ensemble`.

| Function | Returns | Description |
|---|---|---|
| `kf_init_reef(area_m2, init_cover_pct, dhw_out, dhw_cap)` | `Int32` | Initialise reef; write DHW series into caller buffer. Returns `n_timesteps` (75) or `-1`. |
| `kf_run_ensemble(deploy_*, n_runs, covers_out, covers_cap, lower/median/upper_out, stats_cap, n_ts_out, n_valid_out)` | `Int32` | Run ensemble; fill caller buffers in column-major order. Returns `0`, `-1` (uninitialised/error), or `-2` (buffer too small). |

All output buffers are caller-allocated `Float32` arrays. Column-major layout matches Julia's native array order; the Rust side must account for this when reshaping.

**Phase 1 limitation:** `init_cover_pct` and all `deploy_*` parameters are accepted at the ABI boundary but not yet wired into the simulation.

## Usage

Install from the Julia Registry:

```julia
julia> ] add Kora

# Optionally add visualization packages as above
julia> ] add CairoMakie
```

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
run_model!(reef, env)

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
# Fit both growth and survival models together
results = process_ecorrap_models(
    "path/to/ecorrap_adult_juv_combined.csv",
    "path/to/ecorrap_to_species.csv";
    region="offshore_north",
    growth_degree=1,
    survival_degree=2
)

growth_models   = results.growth_fits
survival_models = results.survival_fits

# Or fit them separately
growth_results = process_growth_models(
    "path/to/ecorrap_adult_juv_combined.csv",
    "path/to/ecorrap_to_species.csv";
    region="offshore_north",
    degree=1
)

survival_results = process_survival_models(
    "path/to/ecorrap_adult_juv_combined.csv",
    "path/to/ecorrap_to_species.csv";
    region="offshore_north",
    degree=2
)
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
│   ├── Kora.jl                 # Module entry point
│   ├── stats.jl                # Statistical utilities
│   ├── metrics.jl              # RMSE, R², correlation metrics
│   ├── corals/                 # Coral biology models
│   │   ├── growth_model.jl     # Polynomial growth with habitat constraint
│   │   ├── mortality_model.jl  # DHW/size/depth-dependent bleaching
│   │   ├── recruitment.jl      # Larval production and genetic selection
│   │   └── size_classes.jl     # Size distribution management
│   ├── reefs/                  # Reef-level simulation
│   │   ├── ReefState.jl        # Core state data structure
│   │   ├── reef_dynamics.jl    # Annual timestep dynamics
│   │   ├── run_model.jl        # Single simulation runner
│   │   └── run_ensemble.jl     # Ensemble runner
│   └── interface/              # Data ingestion and model fitting
│       ├── observations.jl     # EcoRRAP data processing
│       ├── regressions.jl      # Growth and survival model fitting
│       └── create_models.jl    # High-level model creation API
├── ext/
│   ├── MakieExt/             # Makie visualization extension
│   └── UnicodePlotsExt/      # Terminal plotting extension
├── bin/
│   └── main.jl               # Interactive web dashboard
├── assets/
│   ├── target_groups.csv     # Functional group definitions
│   └── models/               # Serialized fitted models
└── build/
    ├── bridge_aot.jl         # @ccallable entry points for AOT compilation
    ├── build.ps1             # Build driver for Windows (PowerShell)
    └── build.sh              # Build driver for Linux/macOS
```

## AI Usage Disclosure

AI was used later in the development process to assist in code structure/organization,
documentation setup, and test generation. All AI generated material was reviewed.
