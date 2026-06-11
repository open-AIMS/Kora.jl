# Ensemble Analysis

This tutorial covers running multi-scenario ensembles and interpreting their output. The framing throughout is scenario space exploration: each ensemble member represents one point in a space of uncertain inputs, not a draw from a probability distribution. For the conceptual background, see [Decision Support Under Uncertainty](../concepts/decision-support-under-uncertainty.md).

## Running an Ensemble

`run_ensemble!` accepts a parameter matrix and runs one simulation per column. Each column defines a single scenario with specific ecological and biological assumptions (initial population structure, growth rates, recruitment dynamics). Varying the columns samples different points in the parameter space while holding the environmental scenario constant, which characterises how outcomes shift across assumptions about internal reef dynamics.

The parameter matrix has shape `(n_params, n_members)`. Rows 1 through 16 are population parameters consumed by `set_population!`. If the matrix has more than 16 rows, rows 17 through `16 + n_groups` are per-group location scalers, and the final two rows are `recruits` and `self_seed` values for `run_model!`.

The following example runs 50 scenarios, varying only a DHW scaling factor applied by the population parameter rows.

```julia
using Kora
using Random

rng = Random.default_rng()

reef = initialize_reef(; n_timesteps=50, n_locs=5, area=200.0)
initialize_coral_population!(reef; rng=rng)
environ = generate_example_environment(50, 5; rng=rng)

n_members = 50

# Build a 16-row parameter matrix with one column per scenario.
# Row 1: initial population density, varied across scenarios.
# Rows 2-6: functional group proportions, must sum to 1.0.
# Rows 7-16: size-distribution parameters; zeros here use the default LogNormals.
baseline_params = zeros(Float64, 16, n_members)
baseline_params[1, :] = range(0.5, 2.0; length=n_members)
baseline_params[2:6, :] .= [0.1, 0.2, 0.25, 0.2, 0.25]

results = run_ensemble!(reef, environ, baseline_params; rng=rng)
```

**Note on environmental uncertainty:** This ensemble samples ecological and biological parameter space while holding `environ` constant. To account for deep uncertainty in environmental futures (climate trajectories, thermal adaptation rates), repeat this ensemble procedure with different environmental representations. Each environmental scenario paired with the full parameter matrix produces a complete characterization of outcomes across both ecological and climate uncertainty.

`run_ensemble!` returns a `NamedTuple` with the following fields.

`cover` is an `Array{Float32, 3}` of shape `(n_timesteps, n_locations, n_members)` containing total coral cover in m^2 for each scenario.

`group_cover` is an `Array{Float32, 4}` of shape `(n_timesteps, n_locations, n_groups, n_members)` containing per-functional-group cover.

`juvenile_cover` is an `Array{Float32, 4}` of shape `(n_timesteps, n_locations, n_groups, n_members)` containing sub-mature coral cover.

`wild_dhw_tolerances` is an `Array{Float32, 5}` of shape `(n_timesteps, n_locations, n_groups, 2, n_members)` where the fourth dimension holds the mean (index 1) and standard deviation (index 2) of the population tolerance distribution.

`params` is the input parameter matrix, returned unchanged.

Progress is reported to the terminal during the run.

## Reading Ensemble Output

With the results in hand, per-scenario cover trajectories and summary statistics can be extracted across the scenario space.

```julia
# Total cover across all locations, all timesteps, all scenarios
cover = results.cover   # (n_timesteps, n_locs, n_members)

# Mean cover across locations at each timestep and scenario
using Statistics
mean_cover = dropdims(mean(cover; dims=2); dims=2)   # (n_timesteps, n_members)

# 10th, 50th, 90th percentile ribbons across scenarios at each timestep
# (averaging over locations first)
ts_mean = dropdims(mean(cover; dims=2); dims=2)   # (n_timesteps, n_members)
p10 = [quantile(ts_mean[t, :], 0.10) for t in axes(ts_mean, 1)]
p50 = [quantile(ts_mean[t, :], 0.50) for t in axes(ts_mean, 1)]
p90 = [quantile(ts_mean[t, :], 0.90) for t in axes(ts_mean, 1)]
```

The 10th, 50th, and 90th percentile ribbons summarise scenario-space coverage. They are not a confidence interval. A confidence interval implies a known sampling distribution and makes claims about parameter estimation. These ribbons show where model outcomes cluster across the defined scenarios and where they spread. A narrow ribbon means outcomes are similar across the sampled scenario space. A wide ribbon means the answer is sensitive to which scenario materialises, which is itself useful information for decision support.

## Comparing Two Strategies Across the Ensemble

A core use case is comparing restoration against a no-restoration baseline while holding the climate scenario fixed. This isolates the effect of the intervention.

Both strategies are run with the same parameter matrix. The climate scenario columns are identical across the two runs, so any difference in outcomes is attributable to the deployment schedule rather than to differences in the climate input.

```julia
# Baseline run: no deployments
reef_base = initialize_reef(; n_timesteps=50, n_locs=5, area=200.0)
initialize_coral_population!(reef_base; rng=Random.default_rng())
results_base = run_ensemble!(reef_base, environ, baseline_params; rng=Random.default_rng())

# Restoration run: configure deployment schedule in reef_restore before running
reef_restore = initialize_reef(; n_timesteps=50, n_locs=5, area=200.0)
initialize_coral_population!(reef_restore; rng=Random.default_rng())
# ... set reef_restore.deployment_times here ...
results_restore = run_ensemble!(reef_restore, environ, baseline_params; rng=Random.default_rng())

# Restoration benefit per scenario at the final timestep
base_final   = results_base.cover[end, :, :]     # (n_locs, n_members)
restore_final = results_restore.cover[end, :, :] # (n_locs, n_members)
benefit = restore_final .- base_final             # (n_locs, n_members)
```

The result is the distribution of restoration benefit across the scenario space. Some scenarios will show a large benefit; others will show little. Examining that distribution identifies which types of conditions are required for restoration to make a meaningful difference.

## Introduction to Outcome Partitioning

Outcome Partitioning divides ensemble members into groups based on whether the simulated outcome meets a defined target. The goal is to identify which regions of the scenario space consistently produce acceptable outcomes and which do not.

Start by defining an acceptability criterion. The example below classifies each scenario as acceptable if mean coral cover across all locations exceeds 15% of reef area at year 50.

```julia
area_per_loc = 200.0f0

# Mean cover across locations at final timestep, one value per scenario member
final_cover = dropdims(mean(results.cover[end, :, :]; dims=1); dims=1)

# Express as fraction of reef area
final_cover_frac = final_cover ./ area_per_loc

# Label each scenario: true = acceptable outcome, false = not acceptable
acceptable = final_cover_frac .> 0.15f0
```

With this labelling, the parameter matrix and outcome labels can be exported to CSV for downstream analysis.

```julia
using CSV, DataFrames

df = DataFrame(
    transpose(baseline_params),
    ["param_$i" for i in axes(baseline_params, 1)]
)
df[!, :acceptable] = acceptable

CSV.write("ensemble_labelled.csv", df)
```

With the ensemble labelled in this way, regions of the scenario space that consistently fall into each outcome group can be identified. Tools such as PRIM (Patient Rule Induction Method) or CART can identify the input conditions that best separate the two groups. Kora does not implement these algorithms internally. It provides the ensemble output that feeds into them.

## Choosing Ensemble Size

The right ensemble size depends on how many parameters are varied and how precisely the boundaries between outcome groups must be characterised.

Varying more parameters requires more ensemble members to achieve the same coverage of the scenario space. Characterising outcome group boundaries precisely also requires more members. As a rule, err on the side of more runs rather than fewer.

Kora's annual timestep keeps individual run times short. A 50-scenario ensemble over 5 locations and 75 timesteps completes in seconds on a modern laptop. A 1000-scenario ensemble over the same setup is practical in a single session. That run count is typical for the kind of scenario space exploration that supports Scenario Discovery analysis.

A small ensemble is sufficient to verify setup; run count can be scaled up for final analysis.
