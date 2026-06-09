# Coral Biology Background

> Supporting reading. Not required to run simulations.

This page provides biological context for the modelling choices in Kora. It is intended for readers who want to understand why the model is structured the way it is, not for readers who simply want to run simulations.

## Coral Functional Groups

Kora represents coral communities using five functional groups: tabular Acropora, corymbose Acropora, branching non-Acropora, small massives, and large massives. Each group has distinct growth rates, mortality schedules, and size distributions calibrated from EcoRRAP field data.

Within-group variation in thermal tolerance is large. It reflects both symbiont identity (the algal symbionts in the family Symbiodiniaceae differ substantially in thermal tolerance) and host genotype as independent contributors. A single functional group therefore contains a distribution of tolerance values rather than a single tolerance level. Kora represents this as a truncated-normal distribution across colonies within each group.

## Degree Heating Weeks

Degree Heating Weeks (DHW) is the standard metric for quantifying accumulated thermal stress on coral reefs. It measures the sum of weekly thermal anomalies above the local bleaching threshold, expressed in units of degrees C times weeks. NOAA Coral Reef Watch DHW accumulates only in weeks where sea surface temperature exceeds the monthly mean maximum by at least 1 degree C.

In Kora, DHW values drive bleaching probability and post-bleaching mortality. DHW values are associated with bleaching and mortality as population-level expectations, not deterministic thresholds. At a given DHW level, some fraction of the colony population bleaches and some fraction of those colonies die. The exact fractions depend on the tolerance distribution of the population, which shifts over time through selection.

An important caveat noted by Middlebrook et al. (2008) is that the rate at which thermal stress accumulates affects bleaching outcome independently of the total DHW. A rapid thermal anomaly can cause more bleaching than a slow one of the same total magnitude, because corals have limited acclimation time. Kora uses annual DHW values and cannot resolve within-year accumulation dynamics. This is a known limitation for scenarios involving acute, rapid thermal events.

## Bleaching Biology and the Diameter Abstraction

Bleaching is the temporary expulsion or loss of function of the symbiotic algae that live in coral tissue. When a colony bleaches, its tissue whitens and it loses most of its energy source. Bleached colonies may recover if thermal stress is relieved quickly, or die if stress is prolonged or severe.

In Kora, bleaching outcome is abstracted as a change in colony diameter. A bleaching event reduces the modelled diameter of affected colonies. This is an explicit model simplification. It captures the net effect of the bleaching and partial-recovery process without modelling the recovery phase in detail. The diameter reduction is not empirically calibrated against recovery trajectory data from specific bleaching events. It is a structural assumption about the magnitude of setback that bleaching imposes on colony development.

From a conceptual standpoint, this abstraction could be understood as analogous to super-individuals in agent-based modelling—where a diameter reduction reflects population-level dynamics (mortality, partial recovery, cohort size change) rather than individual coral shrinkage. This framing underscores that the model operates at population scales rather than tracking individual trajectories.

The abstraction is appropriate for the model's purpose. Kora is not designed to predict recovery trajectories at individual colonies. It is designed to capture how bleaching events affect population-level cover dynamics over years to decades.

## Thermal Adaptation and the Selection Assumption

Kora applies the Breeder's equation to advance the tolerance distribution of each functional group after bleaching events. The logic is that bleaching is selective: colonies with lower thermal tolerance bleach and die at higher rates, and survivors' offspring will therefore be more tolerant on average. The shift in mean tolerance per generation depends on the selection differential and the heritability parameter ($h^2$).

The Breeder's equation is a model assumption, not a confirmed GBR outcome. Evidence for heritable thermal tolerance in GBR corals exists -- Csaszar et al. (2010) and Dixon et al. (2015) provide relevant data for GBR populations -- but the degree to which in-situ natural selection is currently advancing thermal tolerance at management-relevant rates remains an active research question. Readers should verify before publication that the Csaszar et al. (2010) citation specifically reports $h^2$ for thermal tolerance traits rather than CO2-stress traits, as the distinction matters.

In Kora, $h^2$ is held constant at $0.3$ across all timesteps, although this could be made a configurable parameter in future releases. Three caveats apply.

First, heritability is population- and environment-specific. A value appropriate for one GBR population may not apply to another, and heritability estimates can change as environmental conditions shift.

Second, the Breeder's equation is strictly a single-generation prediction. Applying it iteratively over many generations erodes additive genetic variance over time. Holding $h^2$ constant across many timesteps is therefore an explicit modelling assumption that overestimates the rate of adaptation in long projections. This is acceptable for the near-term scenario comparisons Kora is designed for, but should be flagged in analyses extending beyond a few decades.

Third, $h^2 = 0.3$ is the default value. It can be varied as a scenario axis. Sensitivity analysis across a range of h^2 values is recommended for any analysis where adaptation rate is a key driver of outcomes.

## Annual Timestep and Temporal Limitations

The annual timestep is appropriate for scenario-space exploration at decadal scales. It is not appropriate for questions about sub-annual dynamics.

The following processes cannot be resolved at annual resolution.

- Within-year bleaching timing and duration.
- The rate at which thermal stress accumulates during a bleaching event.
- Back-to-back event compounding when consecutive bleaching events occur less than one year apart.

The last point is particularly relevant to the recent GBR record. The GBR experienced severe mass bleaching events in 2016, 2017, 2020, and 2022. The 2016 and 2017 events were separated by less than one year. A model with annual resolution cannot represent the compounding mortality that occurs when a reef that has not recovered from one event is hit by another. Users modelling conditions analogous to the 2016-2017 sequence should interpret Kora outputs with this limitation in mind.

## Depth Refugia and Connectivity Assumptions

Kora includes a depth coefficient that reduces bleaching impact on deeper corals. The rationale is that thermal anomalies are attenuated with depth. Cooler deeper water can buffer colonies from surface temperature spikes, at least under moderate warming scenarios.

The depth protection effect decreases as background ocean temperatures rise. As the ocean warms over decadal timescales, the thermal refuge value of deeper water diminishes because deep water itself warms. The depth refuge assumption is most defensible for near-term moderate warming scenarios. It should be applied with caution to late-century high-emission projections, where mesophotic reef thermal refuge capacity may be substantially reduced.

The self-seeding proportion parameter controls what fraction of larvae settle locally on the reef where they were produced. The remainder are assumed to be lost from the simulated system (or to contribute to other reefs not being modelled). The default value is based on general estimates of larval retention at reef scales. It should be validated against site-specific connectivity data where available. Reefs with high isolation or strong local hydrodynamics that promote retention may warrant higher values. Reefs embedded in well-connected reef networks may warrant lower values, with the understanding that external larval supply from neighbouring reefs is not currently modelled.
