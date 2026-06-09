# Decision Support Under Uncertainty

This page provides the intellectual framework behind how Kora is used. The [What Can Kora
Tell Me?](../what-can-kora-tell-me.md) page is a shorter action-oriented introduction.


## Prediction vs. Scenario Analysis

A prediction asks what will happen. It requires a reliable model of the future and assumes that the model can be validated against a known outcome distribution. A scenario analysis asks what would happen if conditions were X. It requires only a plausible forward model that correctly represents the mechanisms connecting inputs to outputs.

Kora is a model for scenario assessment. The annual timestep and functional-group aggregation are not approximations to some higher-fidelity truth. They are the appropriate fidelity for the question being asked. When the question is how restoration outcomes compare across a wide range of climate and ecological conditions, adding sub-annual resolution would increase computational cost and model complexity without improving the quality of the answer.

This distinction matters for how Kora outputs should be interpreted. An ensemble run does not produce a predictive probability distribution that represents what will actually occur. Instead, it produces an empirical map showing how coral cover outcomes vary across a structured scenario space defined by different climate trajectories, ecological parameters, and deployment strategies. The observed distribution of outcomes reflects variation in scenario assumptions, not the likelihood of different futures. The analyst's job is to examine that map and ask which strategies perform acceptably across the range of scenarios considered most plausible or concerning.

## Deep Uncertainty in Coral Reef Management

Not all uncertainty is the same. Some parameters are uncertain but bounded by data. Others are genuinely unknown in a stronger sense: we cannot assign reliable probabilities to their future values even in principle, because the underlying processes are not stable enough or well-enough understood.

The following inputs to reef management decisions fall into the deep uncertainty category.

- Future emission trajectories and their translation into sea surface temperature anomalies at specific reef locations.
- Local thermal anomaly evolution at management-relevant spatial scales, which is not resolved by global climate models.
- The degree of in-situ thermal adaptation on the Great Barrier Reef under ongoing selection pressure.
- Larval connectivity between reefs at the spatial scales relevant to restoration planning.

The following parameters are uncertain but not deeply so.

- Polynomial growth model parameters, which are well constrained by EcoRRAP field data.
- Background mortality rates, which are empirically bounded by long-term monitoring.

Deep uncertainty means we cannot assign reliable probabilities to future states. It is not the same as imprecision in a well-defined distribution. When a parameter is uncertain in the ordinary sense, standard sensitivity analysis and probabilistic framing are appropriate. When a parameter is deeply uncertain, those tools can mislead by implying more structure than exists. The appropriate response is to evaluate strategies across many possible values of that parameter rather than to estimate a best guess.

## Robustness as the Design Target

A robust strategy performs adequately across many futures rather than optimally in one. Consider two strategies: Strategy A produces a 20% cover gain in every scenario examined. Strategy B produces a 40% cover gain in the median scenario but only 5% gain under high-DHW conditions. Strategy A is more robust, and in most management contexts that robustness is more valuable than B's upside.

Kora ensemble output makes this comparison directly. The ensemble does not provide an expected outcome. It provides the distribution of outcomes across the scenario space, which is needed to evaluate robustness.

This framing changes how results should be communicated. The relevant question is not "what will coral cover be in 2050?" but "under which conditions does this deployment strategy maintain cover above X%, and how likely are those conditions relative to the scenarios we are concerned about?" Kora provides the first part of that answer. The second part is a matter of judgment about which scenarios are plausible.

Robustness analysis also changes what counts as a good result. A strategy that dominates under most scenarios but fails badly under a small number of them may be less desirable than one that never fails badly, even if its average performance is lower. This is sometimes called satisficing: finding strategies that clear a threshold across many futures rather than maximising expected performance.

## Scenario Discovery

Scenario Discovery is a method for identifying the specific conditions under which a management strategy succeeds or fails. Rather than asking which strategy is best on average, Scenario Discovery asks what would have to be true, about climate trajectories, parameter values, or deployment scales, for a given strategy to achieve its target outcome.

Kora supports Scenario Discovery by acting as the simulation engine that generates the scenario outcomes. Three features make Kora suitable for this role.

The annual timestep makes large ensembles computationally practical. Running thousands of scenarios is feasible where a daily-timestep model would not be.

The structured parameter space, including DHW scenario, heritability, deployment scale, self-seeding proportion, and others, provides a representative spread of ecological and climate outcomes. Each dimension of that space can be treated as a Scenario Discovery axis.

Outputs from Kora ensemble runs are well-suited to workshop discussions. Managers can explore what would have to be true for a strategy to work, rather than being asked to accept or reject a single projection.

Kora does not implement Scenario Discovery algorithms internally. It provides the ensemble outputs that feed into external Scenario Discovery tools. See the [Ensemble Analysis tutorial](../tutorials/ensemble-analysis.md) for how to run an ensemble and export results for downstream analysis.
