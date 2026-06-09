
# Model Overview {#Model-Overview}

This page describes what Kora simulates and the assumptions behind each model component. It is intended for analysts who have read [What Can Kora Tell Me?](what-can-kora-tell-me.md) and want to understand the model structure before running it. If you are ready to run your first simulation, go to [Getting Started](getting-started.md).

## Simulation Structure {#Simulation-Structure}

Kora uses an annual timestep. Each run covers a configurable number of years, typically 50 to 100. Reef locations are treated as independent units: there is no connectivity between locations in the base model. Every location runs through the same sequence of operations at each timestep.

Within each timestep the following operations are applied in order.
1. Bleaching mortality is applied based on the DHW value for that location and year.
  
2. Background survival is applied, removing colonies stochastically according to size-dependent survival curves.
  
3. Surviving colonies grow according to polynomial growth functions.
  
4. New recruits are added to each location based on larval production from mature colonies.
  

The five functional groups represented in the bundled models are tabular Acropora, corymbose Acropora, branching non-Acropora, small massives, and large massives.

## Coral Representation {#Coral-Representation}

Kora does not track individual colonies. Each functional group at each location is represented as a collection of colony diameters stored in a vector. The distribution of those diameters approximates a truncated log-normal distribution, with a mean and spread that shift over time as bleaching, survival, growth, and recruitment alter the population.

Cover is computed from diameters as the sum of colony projected areas, where each colony contributes $\frac{\pi}{4}\left(\frac{d}{100}\right)^2 \text{ m}^2$ for a diameter $d$ in centimetres.

Bleaching susceptibility differs between functional groups, with Acropora groups generally more susceptible than massive corals. Within-group variation in susceptibility is large and reflects both symbiont identity and host genotype. This variation is represented implicitly through the size-dependent susceptibility function: smaller colonies are more susceptible than larger ones.

## Thermal Stress and Bleaching {#Thermal-Stress-and-Bleaching}

DHW inputs drive bleaching mortality. The bleaching function uses the DHW value for a given location and timestep together with the current population tolerance distribution to compute the proportion of the population affected. Below 4 deg-weeks, no bleaching mortality is applied. Above 4 deg-weeks, the proportion of the population affected increases with DHW. DHW values of 8 or more deg-weeks are associated with significant bleaching and mortality in the model. These are population-level expectations, not deterministic thresholds.

The depth of the reef location modifies the bleaching effect. Deeper locations experience reduced bleaching for a given DHW value, following the relationship described by Baird et al. (2018). Depth is set per location via the `depths` argument to `initialize_reef`.

Bleaching reduces colony diameter, which reduces cover. Colonies below a size threshold experience proportionally greater mortality. The affected fraction of the population is removed or reduced in size; surviving colonies retain their diameters and continue to grow in subsequent timesteps.

## Thermal Adaptation {#Thermal-Adaptation}

The model includes an optional selection mechanism based on a form of the Breeder's equation. If thermal tolerance heritability ($h^2$) is non-zero, each bleaching event shifts the mean thermal tolerance of the surviving population upward. Colonies that survived the bleaching event had higher-than-average tolerance, and that advantage is partially heritable in the next generation.

This mechanism is an explicit modelling assumption. It represents a plausible pathway for in-situ thermal adaptation on the Great Barrier Reef, but the degree to which it occurs at management-relevant rates is not established. Heritability is held constant across timesteps, which is a deliberate simplification. The parameter space for heritability is part of the scenario space that ensemble runs are designed to explore.

## Restoration Deployments {#Restoration-Deployments}

Corals can be deployed at specified timesteps and locations. Deployments add colony diameters directly to the size distribution of the target functional group. Users configure deployment schedules before running the simulation by populating `reef_state.deployment_times`. The model applies deployments internally during the run.

Deployed corals and wild corals are tracked in separate population arrays (`wild_population` and `deployed_population`), which allows their contributions to total cover to be distinguished in the output.

## Growth and Survival Models {#Growth-and-Survival-Models}

Growth and survival are represented as polynomial functions fitted to EcoRRAP survey data from offshore northern Great Barrier Reef sites. There is one growth function and one survival function per functional group. The bundled models are loaded automatically when the package is loaded.

Growth is a function of current colony diameter. The polynomial maps $\log(d)$ to expected diameter increment per year. Growth slows as reef-wide cover approaches the carrying capacity, which is determined by reef area and the `density` parameter.

Survival is modelled using a complementary log-log regression of colony diameter against annual survival probability. Larger colonies have higher survival probability.

Users who need region-specific parameters can supply alternative model sets using `load_models`. The `growth_models` and `survival_models` arguments to `initialize_reef` accept any fitted model collection that implements the required interface.

## Scope and Limitations {#Scope-and-Limitations}

The following processes are not represented in the current model: crown-of-thorns starfish outbreaks, cyclone damage, flood plumes and water quality effects, disease, macroalgae competition, and local anthropogenic stressors. The model is calibrated for scenario-space exploration at annual resolution. It is not designed for sub-annual dynamics, mechanistic ecological research, or sites where non-DHW disturbances dominate reef state.

See [What Can Kora Tell Me?](what-can-kora-tell-me.md) for a summary of the question types the model is and is not suited to answer.
