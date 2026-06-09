# Model Calibration

> This page provides a high-level overview of the calibration workflow. Detailed
> step-by-step documentation, worked examples, and a stable API are still being
> developed.

## What Calibration Does

Kora's default bundled models and parameter ranges are derived from offshore northern
Great Barrier Reef survey data. When Kora is applied to a different reef or a specific
site with its own monitoring record, those defaults may not reproduce the observed
dynamics well. Calibration addresses this by searching for parameter sets that minimise
the mismatch between model output and the observed coral cover time series at the target
site.

The result of a successful calibration is not a single best-fit parameter vector. It is
an *ensemble* of parameter vectors, typically on the order of 100 members, that all
reproduce the observations within an acceptable tolerance. This ensemble captures the
range of parameter combinations consistent with the data and is used as the starting
point for subsequent forward projections and sensitivity analysis.

## What Gets Calibrated

The calibration adjusts 23 parameters that together control the initial reef state and
the demographic behaviour of the five functional groups. They fall into five categories.

| Category | Parameters |
|---|---|
| Colony density | Maximum colony density ($\text{colonies}/\text{m}^2$) |
| Functional group proportions | Relative abundance of each of the five groups at initialisation |
| Initial size distribution | Mean and standard deviation of colony diameter for each group |
| Growth scalers | Per-group multipliers applied to the fitted polynomial growth functions |
| Recruitment and connectivity | External larval supply rate and self-seeding proportion |

Functional group proportions are sampled from a Dirichlet distribution (via a
gamma-to-Dirichlet transform) to ensure they sum to one. All other parameters are
sampled uniformly within bounds derived from the site configuration.

## Calibration Approach

The calibration is formulated as a black-box optimisation problem. An objective
function runs a full Kora simulation for a candidate parameter vector and returns
the root-mean-square error (RMSE) between the simulated total coral cover time series
and the observed cover from LTMP or EcoRRAP monitoring data at the target site. A
lower RMSE indicates a better-fitting candidate.

Multiple independent optimisation trials are run in parallel. Candidates below a
fitness threshold are collected into a shared pool, and the best-performing members
are retained as the calibrated ensemble. Because the objective landscape is
multi-modal, with many different combinations of growth scalers, densities, and
proportions reproducing a similar cover trajectory, so running many trials and
collecting a diverse ensemble is more informative than finding a single global
optimum.

## What the Calibration Produces

After calibration, the following outputs are available.

**Ensemble parameter matrix.** A matrix of shape $(23 \times N)$ where each column
is one calibrated parameter vector and $N$ is the ensemble size (typically 100).

**Performance metrics.** For the ensemble as a whole and for the single best-fit
member: RMSE, Pearson correlation, Kendall tau, bias ($\beta$), and a variability
ratio ($\alpha$) comparing error variance to observed variance.

**Parameter identifiability summary.** Per-parameter spread statistics (coefficient
of variation, relative median absolute deviation) that indicate how tightly the data
constrain each parameter. Poorly identified parameters have a wide posterior range
relative to the prior bounds; well-identified parameters are strongly narrowed.

## Current Status

This workflow is currently implemented as a research script and is not yet available
as a built-in Kora function. The planned direction is to expose it as either a
standard method within Kora or as a separate, model-agnostic calibration tool that
can be applied to any model that accepts a parameter vector and returns a scalar
fitness value.

See [Ensemble Assessment](ensemble-assessment.md) for the next step: analysing the
calibrated ensemble to understand parameter sensitivity and uncertainty.
