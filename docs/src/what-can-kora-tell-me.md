# What Can Kora Tell Me?

Kora is a simulation model for decision support in coral reef restoration planning. It is built for conditions of deep uncertainty, where the future trajectory of climate stress, larval supply, and thermal adaptation cannot be reliably predicted, and the goal is to find strategies that perform acceptably across a wide range of futures rather than to optimise for a single expected outcome.

## Decision Context Kora Supports

Reef managers and ecologists typically need to make decisions before the relevant science is resolved. Kora is designed for that setting. It is not a forecasting tool. It is a tool for comparing strategies under conditions where the future cannot be reliably predicted.

These are questions about robustness under uncertainty, not questions asking for a single-number prediction.

The following are representative questions that Kora is designed to help answer.

- "Should we deploy heat-tolerant corals at this location, and if so, how many colonies and over what reef area?"
- "How do different climate change scenarios affect potential restoration outcomes?"
- "How sensitive is the outcome to our assumptions around thermal adaptation and coral demographic rates?"

The shift in question type matters. Asking which strategy is best requires predicting the future. Asking which strategy avoids poor outcomes across all the climate scenarios we are uncertain about does not. Kora is intended for the second type of question.

## Required Inputs

Kora requires a small set of inputs to run. The minimum set consists of DHW time series provided at annual resolution, reef area and an initial colony size distribution, and a deployment schedule if a restoration intervention is being modelled. Group-specific biological parameters are optional and have defaults derived from EcoRRAP field data.

Real DHW data can be wrapped for use in Kora using the `generate_environment()` function. See the Input Data Reference for full format details.

## What Outputs Look Like

Kora produces time series of coral cover, colony count, and size-class distributions for each location and functional group. Bleaching event severity is recorded per timestep. When the model is run as an ensemble, summary ribbons at the 10th, 50th, and 90th percentile are the standard output format. Restoration and baseline runs can be compared directly.

Outputs are best understood as inputs to a decision, not as predictions of a known future. Each ensemble run produces a separate time series. The built-in `viz.ensemble_timeseries` plot summarises these as a median line with a quantile envelope. Where the envelope is narrow, outcomes are similar across scenarios. Where it is wide, the answer depends strongly on which climate or parameter scenario materialises, and that spread is itself decision-relevant information.

## Is This Tool Right for My Question?

The following categories describe where Kora is well-suited, where it can be used with care, and where it is less appropriate.

Kora is well-suited for the following uses.

- Comparing restoration strategies under a range of climate scenarios.
- Assessing sensitivity of outcomes to deployment scale and timing.
- Identifying conditions under which a strategy succeeds or fails, sometimes called Scenario Discovery.
- Tolerance-uplift cost-benefit framing under deep uncertainty.

The following uses require care.

- Short-term site-specific projections: the model can provide an indication provided no major non-DHW disturbances are expected at the site during the forecast period.
- Mechanistic ecological research or validating ecological theory: Kora is designed for decision support at management scales, not for resolving fine-grained ecological mechanisms.

Kora is less suited for the following situations.

- Systems where crown-of-thorns, cyclone, or water quality dynamics are the dominant driver of reef state at the time of interest.

Crown-of-thorns, cyclones, flood plumes, disease, macroalgae competition, and local stressors are not currently represented. The model's growth and survival components are structured around a dispatch system that can be extended per functional group. Disturbance processes are currently implemented on a case-by-case basis. This is a known extensibility target rather than a fundamental constraint.

## Appropriate Precision

The annual timestep and functional-group aggregation in Kora are deliberate design choices, not limitations to be apologised for. The question Kora is designed to answer is about the direction and magnitude of restoration benefit across many futures. Resolving sub-annual dynamics at a single site is a different question that calls for a different model.

The EcoRRAP empirical grounding is the basis for trusting Kora's outputs for their intended purpose. The model parameters are constrained by field observations rather than assumed from first principles. This means the model is well-calibrated for the scale and question type it was designed for, even though it is not a high-fidelity mechanistic simulator.

See the Input Data Reference for guidance on how to provide real DHW data and site-specific parameters.

## How Kora and ADRIA Can Work Together

Kora and ADRIA can be positioned as complementary components in one decision-support stack.

- Kora provides the coral ecology representation. It maps ecological trajectories under alternative climate, adaptation, and restoration assumptions.
- ADRIA provides the downstream assessment layer for scenario ensembles, sensitivity analysis, and comparison of strategy performance across plausible futures.

For non-modellers, this can be read as: Kora simulates coral outcomes under different assumptions, and ADRIA helps compare those outcomes when selecting strategies.

Current integration is workflow-level. Kora outputs can be used in ADRIA workflows, but there is not yet a unified package-level interface. An open Kora.jl issue proposes adopting data structures from open-AIMS/ADRIAIndicators.jl to improve interoperability. This indicates planned direction rather than completed integration.

A practical integration pattern is:

1. Define ecological assumptions and intervention strategies in Kora.
2. Run scenario ensembles in Kora and export outcome metrics such as coral cover trajectories, group-specific outcomes, and restoration-versus-baseline differences.
3. Use ADRIA to compare alternatives against decision objectives, assess sensitivity to assumptions, and evaluate robustness across plausible futures.
4. Update Kora scenarios based on ADRIA findings and repeat until strategy choices are decision-ready.

This separation of roles keeps ecological process representation explicit in Kora while using ADRIA as the analysis layer for robustness and trade-off assessment. Conclusions remain conditional on the modeled processes and scenario space.

See the [ADRIA.jl repository](https://github.com/open-AIMS/ADRIA.jl) and [ADRIAIndicators.jl](https://github.com/open-AIMS/ADRIAIndicators.jl).
