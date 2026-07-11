# Kora and ADRIA in a Decision Workflow

Kora and ADRIA serve different roles in a shared decision-support workflow.

Kora provides the coral ecology representation. It simulates coral population dynamics and resulting cover trajectories under alternative climate, thermal adaptation, and restoration assumptions. These outputs are scenario-conditioned ecological outcomes, not forecasts of the most likely reef future.

ADRIA.jl provides the surrounding assessment layer for scenario ensembles, sensitivity analysis, and comparison of strategy performance across plausible futures. In this context, robustness means identifying options that remain acceptable, or fail less severely, across a defined set of plausible futures.

[CoralBlox.jl](https://github.com/open-AIMS/CoralBlox.jl) is ADRIA's companion coral model. ADRIA is designed to be model-agnostic with respect to ecological model outputs, so Kora can be substituted in any workflow where its representation of thermal exposure, adaptation, and restoration is better suited to the question. One practical distinction: environmental drivers such as DHW are provided to CoralBlox by ADRIA, whereas Kora handles environmental forcing internally.

Current integration is workflow-level: Kora outputs can be transferred into ADRIA analyses, but there is not yet a unified package-level interface. An open Kora.jl issue proposes adopting data structures from open-AIMS/ADRIAIndicators.jl to reduce translation effort and improve interoperability. That alignment is planned work, not a completed interface.

Kora defines ecological scenarios and simulates coral outcome trajectories.
ADRIA structures scenario ensembles, assesses sensitivity to assumptions, and compares decision performance across plausible futures.
Insights are fed back into Kora scenario design for iterative refinement of ecological assumptions and intervention options.

See the [ADRIA.jl repository](https://github.com/open-AIMS/ADRIA.jl), [CoralBlox.jl](https://github.com/open-AIMS/CoralBlox.jl), and [ADRIAIndicators.jl](https://github.com/open-AIMS/ADRIAIndicators.jl).
