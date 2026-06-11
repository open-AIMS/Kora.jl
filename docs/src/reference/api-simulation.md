# Simulation API

```@docs
run_model
run_model!
run_ensemble!
set_population!
assign_scalers!
```

## Deployment (Internal)

`deploy_corals!` is called internally by `run_model!` when coral outplanting is
active. It is documented here for completeness; deployment should be configured through the `run_model!` interface rather than by calling `deploy_corals!` directly.

```@docs
deploy_corals!
```
