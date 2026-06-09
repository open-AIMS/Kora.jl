# Ensemble Assessment

> This page provides a high-level overview of post-calibration ensemble assessment.
> Detailed step-by-step documentation, worked examples, and a stable API are still
> being developed.

## Purpose

Once a calibrated ensemble of parameter vectors is available (see
[Model Calibration](model-calibration.md)), the next question is: *which of those
23 parameters actually matter, and when do they matter?* Ensemble assessment answers
this through a combination of global and constrained sensitivity analysis, temporal
analysis, and uncertainty quantification.

The outputs help identify which parameters are driving outcomes during critical periods
such as bleaching events and post-disturbance recovery, and which parameters the
calibration data has successfully constrained versus left largely unconstrained.

## Global versus Constrained Sensitivity

For comparison, two sensitivity analyses are run on the calibrated ensemble.

**Global (unconstrained) sensitivity** evaluates parameter importance over the full
prior parameter space defined by the original search bounds. It answers: *with no prior
system knowledge and parameters sampled at random from the prior, which
parameters could have the largest effect on coral cover?*

**Constrained (posterior) sensitivity** repeats the analysis using only the parameter
ranges spanned by the calibrated ensemble. It answers: *within the part of parameter
space that fits the observations, which parameters are still driving variation in
outcomes?*

Comparing the two analyses reveals which parameters were successfully constrained by
the calibration (their posterior sensitivity is lower than their prior sensitivity)
and which remain influential even after fitting (they are candidates for further data
collection or scenario axes in forward projections).

Both analyses use the PAWN sensitivity method, which is a distribution-based
global sensitivity index that does not require a specific model structure and makes
no assumptions about linearity or parameter independence.

## Temporal Sensitivity

A single aggregate sensitivity index collapses the time dimension and can miss
important structure. Bleaching events create sharp transitions in which parameters
drive cover dynamics, and the dominant parameters during an acute bleaching year are
different from those that drive recovery in the years that follow.

Temporal sensitivity analysis computes PAWN indices at each timestep independently,
producing a matrix of sensitivity values over parameters and time. This shows, for
example, that growth scalers become more important in recovery years while initial
density and proportion parameters dominate the early simulation period.

A related analysis computes *lagged* sensitivity: for each lag $k \in \{1, 2, 3, 5\}$
years, PAWN is computed on cover at time $t + k$ as a function of parameters at time
$t$. This identifies which parameters drive recovery at different time horizons after
a disturbance.

## Per-Functional-Group Sensitivity

All of the above analyses can be repeated separately for each of the five functional
groups rather than for total coral cover. This reveals group-specific dynamics: the
parameters that control tabular Acropora recovery after bleaching may differ from
those controlling large massive persistence, because the two groups have distinct
growth and mortality schedules.

## Uncertainty Quantification

For each parameter, two quantities are compared.

**Prior standard deviation** is the spread of values sampled uniformly across the
original search bounds.

**Posterior standard deviation** is the spread of values in the calibrated ensemble.

The ratio of these quantities, expressed as a percentage reduction, measures how much
the calibration data narrowed uncertainty for each parameter. Parameters with large
reductions were effectively constrained by the observations. Parameters with small
reductions remain uncertain and are candidates for sensitivity analysis axes in
forward projections.

## Current Status

This workflow is currently implemented as a research script and is not yet available
as a built-in Kora function. The planned direction mirrors that of the calibration
step: to provide these analyses as either standard Kora methods or as part of a
separate, model-agnostic post-calibration toolkit.
