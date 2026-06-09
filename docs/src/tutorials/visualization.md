# Visualizing Simulation Results

This tutorial covers Kora's built-in visualization functions powered by Makie. These functions create publication-quality figures and can be displayed interactively or saved to disk.

## Basic Timeseries Plot

For a single simulation run, `timeseries()` creates a comprehensive plot showing total coral cover and DHW conditions:

```julia
using Kora
using Random

reef = initialize_reef(; n_timesteps=50, n_locs=10, area=200.0, density=15)
initialize_coral_population!(reef; rng=Random.default_rng())
environ = generate_example_environment(50, 10)

run_model!(reef, environ; rng=Random.default_rng())

fig = Kora.viz.timeseries(reef, environ)
display(fig)
```

The figure shows total coral cover in the top panel and DHW at a specified location in the bottom panel. By default it displays location 1; pass `loc=N` to plot a different location.

To save the figure:

```julia
save("coral_timeseries.png", fig)
```

## Ensemble Timeseries

When running multiple scenarios, `ensemble_timeseries()` visualizes the full ensemble with quantile bands showing the spread of outcomes:

```julia
n_scenarios = 50

# Build parameter matrix: vary density across scenarios
params = zeros(Float64, 23, n_scenarios)
params[1, :] = range(10.0, 25.0; length=n_scenarios)

results = run_ensemble!(reef, environ, params; rng=Random.default_rng())

fig = Kora.viz.ensemble_timeseries(reef, results, environ)
display(fig)
```

The output shows:
- Median coral cover (solid line)
- 95% credible interval (shaded band, 2.5th to 97.5th percentile)
- Individual ensemble member trajectories (faint lines)
- DHW conditions (bottom panel)

Customize quantiles by passing `quantiles=[0.1, 0.5, 0.9]` to show 10th, 50th, and 90th percentiles instead.

## Group-Level Cover

To examine how different functional groups respond, extract per-group coverage:

```julia
group_cov = Kora.group_cover(reef)  # Returns (n_timesteps, n_groups) matrix

# Access cover for a specific group at each timestep
tabular_acro_cover = group_cov[:, 1]
```

Then plot manually or use Makie directly for custom visualizations.

## Juvenile Cover

Colonies below the maturity threshold for their group are classified as juveniles. Track juvenile recruitment and development:

```julia
# Get juvenile cover for all groups at timestep 20
jc = Kora.juvenile_cover(reef, 20)  # Returns vector of length n_groups

# Examine juvenile-to-adult ratio
total_cov = Kora.coral_cover(reef, 20)
juv_fraction = sum(jc) / total_cov
```

## Population Size Distribution Animation

Animate the size distribution of colonies over time to see how populations change during the simulation:

```julia
Kora.viz.animate_population(
    reef, environ,
    loc=1,           # Location index
    grp=1,           # Functional group (1 = tabular Acropora)
    filename="population_evolution.gif"
)
```

This creates an animated GIF showing:
- Histogram of colony diameters (top panel)
- Size distribution changes as cohorts grow or die (middle)
- DHW conditions over time (bottom), with a vertical line marking the current timestep

## Multiple Runs and Comparison

To compare outcomes across different scenarios or strategies:

```julia
# Scenario A: baseline
run_model!(reef, environ; rng=Random.default_rng())
cover_a = Kora.coral_cover(reef)

# Scenario B: different DHW
Kora.reset!(reef)
initialize_coral_population!(reef; rng=Random.default_rng())
environ_alt = generate_example_environment(50, 10)  # Different synthetic DHW
run_model!(reef, environ_alt; rng=Random.default_rng())
cover_b = Kora.coral_cover(reef)

# Plot both timeseries on the same axes
using CairoMakie
fig, ax = CairoMakie.scatter(cover_a; label="Baseline")
lines!(ax, cover_b; label="High DHW", color=:red)
axislegend()
display(fig)
```

## Performance Model Plots

When fitting custom growth or survival models, `growth_performance_plots()` and `survival_performance_plots()` visualize model accuracy against empirical data:

```julia
growth_fits = Kora.fit_growth_models(growth_groupings; degree=2)

fig = Kora.viz.growth_performance_plots(
    growth_groupings,
    growth_fits;
    save_path="growth_model_fits.png"
)
```

This shows fitted curves overlaid on data bins for each functional group, with performance metrics (RMSE, R-squared, correlations) displayed.

## Saving and Exporting

All visualization functions support flexible output:

```julia
# Display interactively (default)
display(fig)

# Save as PNG (default format)
save("output.png", fig)

# Save as SVG (vector format for publication)
save("output.svg", fig)

# Save as PDF
save("output.pdf", fig)

# Specify DPI for raster formats
save("output_highres.png", fig; px_per_unit=2)
```

## Notes on Visualization Limits

The current Makie-based visualization supports plotting single locations or aggregated results. For detailed spatial analysis across many locations, consider exporting raw data and visualizing with external tools (e.g., mapping libraries for geographic display).

Bootstrap confidence intervals are also available via `bootstrap_cover()`, `bootstrap_ensemble_timeseries()` and related functions, but are typically used for internal analysis rather than publication figures.
