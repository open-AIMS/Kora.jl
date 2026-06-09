"""
Standalone script: generates static PNG/GIF figures for the documentation.
Run once from the repo root:

    julia --project=docs docs/generate_figs.jl

Outputs are written to docs/src/assets/.
"""

using Random
using CairoMakie
using Kora

const OUT = joinpath(@__DIR__, "src", "assets")
mkpath(OUT)

rng = Random.default_rng()

# ---------------------------------------------------------------------------
# Figure 1: timeseries (single run)
# ---------------------------------------------------------------------------
reef = initialize_reef(; n_timesteps=50, n_locs=1, area=300.0, density=9)
initialize_coral_population!(reef; rng)
environ = generate_example_environment(50, 1; rng)
run_model!(reef, environ; rng)

fig1 = Kora.viz.timeseries(reef, environ)
save(joinpath(OUT, "viz_timeseries.png"), fig1)
@info "Saved viz_timeseries.png"

# ---------------------------------------------------------------------------
# Figure 2: ensemble_timeseries
# ---------------------------------------------------------------------------
n_scenarios = 50

# Row 1: initial population density (colonies/m²), varied across scenarios.
# Rows 2-6: functional group proportions, must sum to 1.0.
# Rows 7-16: size-distribution parameters (left at 0 -> default LogNormals used).
const DEFAULT_PROPORTIONS = [0.2, 0.2, 0.2, 0.2, 0.2]
params = zeros(Float64, 6, n_scenarios)
params[1, :] = range(5.0, 20.0; length=n_scenarios)
params[2:6, :] .= DEFAULT_PROPORTIONS

results = run_ensemble!(reef, environ, params; rng)

fig2 = Kora.viz.ensemble_timeseries(reef, results, environ)
save(joinpath(OUT, "viz_ensemble_timeseries.png"), fig2)
@info "Saved viz_ensemble_timeseries.png"

# ---------------------------------------------------------------------------
# Figure 3: population animation (single representative frame saved as PNG)
# ---------------------------------------------------------------------------
Kora.viz.animate_population(
    reef, environ,
    1, 1;
    filename=joinpath(OUT, "viz_population_animation.gif")
)
@info "Saved viz_population_animation.gif"

@info "All figures written to $OUT"
