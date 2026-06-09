# Restoration Scenarios

## Overview

This tutorial shows how to use Kora to compare a baseline trajectory against a restoration
intervention across a range of climate scenarios. The goal is not to predict what will happen
at a specific site. The goal is to understand the distribution of restoration benefit across
the climate futures we cannot rule out.

See [What Can Kora Tell Me?](../what-can-kora-tell-me.md) for a plain-language framing of
the questions Kora is designed to answer. The focus here is on the practical steps for
setting up and interpreting restoration comparisons.

The central idea in this tutorial is robustness rather than prediction. A deployment strategy
that consistently produces moderate improvement across many scenarios is often more useful to
a reef manager than one that produces a large improvement under optimistic conditions and
fails under pessimistic ones. Running both baseline and restoration configurations over the
same set of DHW scenarios gives you the distribution of benefit, not a single expected number.

## What a Restoration Configuration Looks Like

In Kora, a restoration scenario is defined by populating the `deployment_times` field of a
`ReefState` object. This field is a three-dimensional array with shape
`(n_timesteps, n_locations, n_groups)`. Each entry holds the number of corals to be deployed
at a given timestep, location, and functional group.

You do not call `deploy_corals!` yourself. It is an internal function executed by the
simulation loop at each scheduled deployment timestep.

The baseline scenario uses the same `ReefState` configuration with `deployment_times` left at
its default value of zero throughout.

**Note:** Full integration of deployment scheduling into the simulation loop is being
finalised in the current development version. The `deployment_times` field and the
`deploy_corals!` function are present in the codebase. Check the package CHANGELOG or
release notes to confirm deployment activation status before relying on this feature in
production workflows. The patterns shown below reflect the intended API.

## Setting Up the Baseline

The baseline represents natural reef dynamics with no intervention. The setup is the
standard three-step sequence of initialisation, population seeding, and running.

```julia
using Random, Kora

n_ts   = 50
n_locs = 10
reef_area = 90.0

environ = Kora.generate_environment(
    my_dhw_matrix;       # Matrix{Float32} of shape (n_ts, n_locs)
)

baseline_reef = Kora.initialize_reef(;
    n_timesteps = n_ts,
    n_locs      = n_locs,
    area        = reef_area,
    density     = 15,
    depths      = 9.0
)

Kora.initialize_coral_population!(baseline_reef; rng = Random.default_rng())
Kora.run_model!(baseline_reef, environ; rng = Random.default_rng())
```

Using the same default RNG state for both the baseline and the restoration run is important.
It ensures that any difference in outcomes between the two runs is attributable to the
deployment schedule rather than to variation in the stochastic processes within the
simulation.

## Setting Up the Restoration Configuration

The restoration configuration starts from an identical initial state. The `deployment_times`
array is filled to describe the intervention schedule before the simulation is run.

```julia
restoration_reef = Kora.initialize_reef(;
    n_timesteps = n_ts,
    n_locs      = n_locs,
    area        = reef_area,
    density     = 15,
    depths      = 9.0
)

Kora.initialize_coral_population!(restoration_reef; rng = Random.default_rng())

# Schedule 50 tabular Acropora deployments at each location in years 5, 10, and 15.
# Group index 1 corresponds to acro_table in the default TARGET_GROUPS ordering.
for deploy_year in [5, 10, 15]
    restoration_reef.deployment_times[deploy_year, :, 1] .= 50.0f0
end

Kora.run_model!(restoration_reef, environ; rng = Random.default_rng())
```

Both runs use the same `environ` object, the same RNG state, and the same initial population.
The only difference is the presence of deployments in the restoration configuration.

## Comparing Outcomes for One Climate Scenario

After both runs complete, compare total coral cover across time.

```julia
# coral_cover returns a Vector{Float32} of length n_timesteps,
# with each element being the total cover summed across all locations at that timestep.
baseline_cover    = Kora.coral_cover(baseline_reef)
restoration_cover = Kora.coral_cover(restoration_reef)

restoration_benefit = restoration_cover .- baseline_cover
```

A positive value at any timestep means the restoration run produced more total cover than
the baseline under that climate scenario. A negative value would indicate the intervention
was counterproductive at that point, which can happen if deploying corals increases
competition for space in a period of low DHW stress.

Plotting `restoration_benefit` as a time series gives a compact summary of when and by how
much the intervention diverges from the baseline under a single climate trajectory.

## Comparing Across Multiple Climate Scenarios

Running both configurations over a set of DHW scenarios rather than a single trajectory
produces the distribution of restoration benefit. This is the step that turns a
single-scenario comparison into a robustness assessment.

```julia
dhw_scenarios = [scenario_1, scenario_2, scenario_3]  # Vector of Matrix{Float32}

cover_baselines    = []
cover_restorations = []

for dhw in dhw_scenarios
    env_i = Kora.generate_environment(dhw)

    b = Kora.initialize_reef(; n_timesteps=n_ts, n_locs=n_locs, area=reef_area, density=15)
    Kora.initialize_coral_population!(b; rng = Random.default_rng())
    Kora.run_model!(b, env_i; rng = Random.default_rng())
    push!(cover_baselines, Kora.coral_cover(b))

    r = Kora.initialize_reef(; n_timesteps=n_ts, n_locs=n_locs, area=reef_area, density=15)
    Kora.initialize_coral_population!(r; rng = Random.default_rng())
    for deploy_year in [5, 10, 15]
        r.deployment_times[deploy_year, :, 1] .= 50.0f0
    end
    Kora.run_model!(r, env_i; rng = Random.default_rng())
    push!(cover_restorations, Kora.coral_cover(r))
end
```

The collection of `cover_baselines` and `cover_restorations` across all scenarios gives you
the spread of outcomes. Each element is a `Vector{Float32}` of total cover over time.
Summarising this spread at representative timesteps, for example year 20 and year 50,
produces the distribution of restoration benefit across the scenarios you considered.

A strategy that produces positive benefit across most scenarios is more robust than one
that only produces benefit in the most optimistic scenario. The width of the benefit
distribution is itself informative: a narrow distribution means the strategy's effectiveness
is not very sensitive to which climate future materialises.

## When Does Restoration Make a Difference?

After running many scenarios, the natural next question is which conditions separate the
scenarios where restoration helped from those where it did not. This is the domain of
Scenario Discovery: instead of summarising average performance, you identify the axes of
the scenario space that predict whether a strategy crosses a defined success threshold.

Inputs to Scenario Discovery might include the DHW trajectory (total accumulated stress over
the simulation period), the initial population state, the deployment scale, and the frequency
of deployment events. Any parameter that varies across your ensemble members is a candidate
axis.

Kora does not implement Scenario Discovery algorithms. It generates the ensemble outputs
that feed into external Scenario Discovery tools. See
[Ensemble Analysis](ensemble-analysis.md) for the full workflow, including how to run a
large ensemble with `run_ensemble!` and export the results for downstream analysis.
