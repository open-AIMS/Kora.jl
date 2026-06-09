# Kora.jl

Kora is a simulation model for decision support in coral reef restoration planning.
It is designed for use under deep uncertainty: it maps how outcomes vary across climate
and ecological scenarios rather than forecasting a single future.

Crown-of-thorns, cyclones, flood plumes, disease, macroalgae competition, and local
stressors are not currently represented. Growth and survival models are extensible per
functional group; disturbances are currently implemented case-by-case and represent a
known extensibility target.

See [What Can Kora Tell Me?](what-can-kora-tell-me.md) for the full scope discussion.

## Quick Start

```julia
using Kora

reef = initialize_reef(n_timesteps=75, n_locs=1, area=60.0, density=10, depths=7.0)
environ  = generate_example_environment(75, 1)
initialize_coral_population!(reef, 1, ceil(Int64, 3 * 60.0); group_proportions=fill(0.2f0, 5))
run_model!(reef, environ)
cover = coral_cover(reef)
```
