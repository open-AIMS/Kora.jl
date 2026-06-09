
# Input Data Reference {#Input-Data-Reference}

This page answers the question: what exactly does Kora expect me to provide, and in what format? It covers the environmental forcing array, the EcoRRAP survey CSV files used for model fitting, the fitted model JSON files, and the full parameter set for `initialize_reef`.

## Section 1: Environmental Forcing (DHW Array) {#Section-1:-Environmental-Forcing-DHW-Array}

The primary environmental input to `run_model!` and `run_ensemble!` is a `YAXArray` produced by `generate_environment` or `generate_example_environment`.

### Structure {#Structure}

|      Property |                                                                                          Value |
| -------------:| ----------------------------------------------------------------------------------------------:|
|          Type |                                                                            `YAXArray{Float32}` |
|   Dimension 1 |                                                       `:timestep` – integer index 1 to n_years |
|   Dimension 2 |                        `:location` – integer index 1 to n_locs; must match `n_locations(reef)` |
|      Variable |                                     `:dhw` – annual Degree Heating Weeks, non-negative Float32 |
|         Units |                              deg C times weeks (standard NOAA Coral Reef Watch DHW definition) |
| Typical range | 0 to 25; values above approximately 20 are at the high end of mid-century SSP5-8.5 projections |


### Construction {#Construction}

To wrap your own DHW data, pass a `Matrix{Float32}` of shape `(n_years, n_locs)` to `generate_environment`.

```julia
dhw_matrix = Float32.(my_dhw_data)   # must be Matrix{Float32}, shape (n_years, n_locs)
environ = Kora.generate_environment(dhw_matrix)
```


To generate synthetic environmental data for testing or exploration, use `generate_example_environment`.

```julia
environ = Kora.generate_example_environment(50, 10)   # 50 years, 10 locations
```


### Input validation and warnings {#Input-validation-and-warnings}

`generate_environment` validates the input matrix and emits warnings for two conditions.

The first warning is triggered when the maximum value in the matrix exceeds 40 deg-weeks. Values in that range are roughly twice the approximately 20 deg-weeks projected under SSP5-8.5 in mid-century scenarios. The warning text indicates this may reflect a data quality issue or an intentionally extreme scenario.

The second warning is triggered when the minimum value across all locations and timesteps exceeds 20 deg-weeks. Real DHW data contains many near-zero values during non-bleaching years and early simulation periods. A floor above 20 deg-weeks across the entire matrix suggests that raw sea-surface temperature (typically 25 to 32 degrees C) may have been passed instead of accumulated degree heating weeks.

Neither warning stops execution. Both are advisory.

### Alignment requirement {#Alignment-requirement}

The number of timesteps and the number of locations in the environment array must exactly match the values used to initialise the corresponding `ReefState`. Passing an environment with a different number of locations than `n_locations(reef)` will produce an error at the start of `run_model!`.

## Section 2: EcoRRAP Survey CSVs {#Section-2:-EcoRRAP-Survey-CSVs}

Two CSV files are required when fitting region-specific growth and survival models. See the [Fitting Models from EcoRRAP Data](../tutorials/fitting-from-ecorrap.md) tutorial for the full fitting workflow.

### Coral demographic survey data {#Coral-demographic-survey-data}

Each row records a single colony observed at two points in time.

|         Column |    Type |                                                 Description |
| --------------:| -------:| -----------------------------------------------------------:|
| `area_t1_sqcm` | numeric |          Colony area at first observation, in $\text{cm}^2$ |
| `area_t2_sqcm` | numeric |         Colony area at second observation, in $\text{cm}^2$ |
|        `taxon` |  string |                                       Species or taxon code |
|     `survival` | integer |                 1 if alive at second observation, 0 if dead |
|   `growth_use` |  string |   "yes" if this row should be used for growth model fitting |
| `survival_use` |  string | "yes" if this row should be used for survival model fitting |
|   `days_t1.t2` | numeric |                           Days between the two observations |
|      `cluster` |  string |            Regional cluster label, e.g. "Offshore_Northern" |
|    `site_code` |  string |       Site identifier used for spatial train/test splitting |


The `cluster` column must be present. Values such as "Offshore_Northern", "Offshore_Central", and "Offshore_Southern" are normalised automatically by `standardize_ecorrap_data!`to`"offshore_north"`,`"offshore_central"`, and`"offshore_south"` respectively.

Column name normalisation is applied automatically. The function `standardize_ecorrap_data!` lowercases all column names and renames `area_t1_sqcm` to `size`, `area_t2_sqcm` to `sizenext`, `taxon` to `taxa`, and `survival` to `surv`. Downstream processing functions expect these normalised names. You do not need to rename columns manually before calling `process_ecorrap_models`, `process_growth_models`, or `process_survival_models`.

### Species-to-functional-group mapping {#Species-to-functional-group-mapping}

Each row maps one taxon code to one functional group.

|         Column |   Type |                                                                         Description |
| --------------:| ------:| -----------------------------------------------------------------------------------:|
|         `Code` | string |                           Taxon code matching the `taxon` column in the survey data |
| `Cscape_group` | string | Functional group identifier; must be one of the five values in `Kora.TARGET_GROUPS` |


The five valid `Cscape_group` values are: `acro_table`, `acro_corym`, `corym_non_acro`, `small_massive`, and `large_massive`. Taxon codes that map to any other string are excluded from model fitting. Taxon codes with no entry in this file are also excluded.

## Section 3: Fitted Model Files (JSON) {#Section-3:-Fitted-Model-Files-JSON}

Growth and survival models are serialised to JSON by `save_models` and deserialised by `load_models`. The same format is used for both the bundled assets and for custom models you fit yourself.

### Bundled model files {#Bundled-model-files}

Kora includes two bundled model files in the package assets directory.

|                                                File |                                             Contents |
| ---------------------------------------------------:| ----------------------------------------------------:|
|   `assets/models/offshore_north_growth_models.json` |   Growth models fitted to offshore northern GBR data |
| `assets/models/offshore_north_survival_models.json` | Survival models fitted to offshore northern GBR data |


These files are loaded automatically when the package starts up. After loading, they are accessible as `Kora.growth_models` and `Kora.survival_models`. Any call to `initialize_reef` that does not supply explicit model arguments uses these defaults.

### File structure {#File-structure}

|            Field |    Type |                                                                                      Description |
| ----------------:| -------:| ------------------------------------------------------------------------------------------------:|
| `format_version` | integer |                                                                      Schema version; currently 1 |
|     `model_kind` |  string |                                                                    Either "growth" or "survival" |
|         `region` |  string |                                                             Optional label supplied at save time |
|      `fitted_at` |  string |                                         ISO 8601 timestamp recording when the models were fitted |
|         `models` |   array | One entry per functional group, each containing polynomial coefficients and numeric range limits |
|    `performance` |  object |                Train and test metrics per group: RMSE, R-squared, Pearson, Spearman, and Kendall |


The `format_version` field is checked on load. An informative error is raised if the version in the file is not in the list of versions supported by the installed Kora release.

### Loading custom model files {#Loading-custom-model-files}

```julia
using Kora

gm = Kora.load_models("my_region_growth.json")
sm = Kora.load_models("my_region_survival.json")

reef = Kora.initialize_reef(;
    n_timesteps     = 50,
    n_locs          = 20,
    area            = 90.0,
    growth_models   = gm,
    survival_models = sm
)
```


### Verifying that a model pair came from the same fitting run {#Verifying-that-a-model-pair-came-from-the-same-fitting-run}

```julia
Kora.check_model_pair_skew("my_region_growth.json", "my_region_survival.json")
```


If the `fitted_at` timestamps in the two files differ by more than 86400 seconds (24 hours), a warning is emitted. This check is advisory: it helps detect accidental mixing of model files from separate fitting runs on different data subsets.

## Section 4: initialize_reef Parameter Reference {#Section-4:-initialize_reef-Parameter-Reference}

`initialize_reef` allocates a `ReefState` sized for the simulation. All population arrays are empty after this call. Call `initialize_coral_population!` before running the simulation.

```julia
reef = Kora.initialize_reef(;
    n_timesteps     = 75,
    n_locs          = 100,
    group_names     = Kora.TARGET_GROUPS,
    density         = 20,
    area            = 90.0,
    depths          = 9.0,
    growth_models   = Kora.growth_models,
    survival_models = Kora.survival_models
)
```


### Parameters {#Parameters}

|         Parameter |                              Type |                Default |                                                                                                                  Description |
| -----------------:| ---------------------------------:| ----------------------:| ----------------------------------------------------------------------------------------------------------------------------:|
|     `n_timesteps` |                             `Int` |                     75 |                                                                                      Number of annual time steps to allocate |
|          `n_locs` |                             `Int` |                    100 |                                                                                                     Number of reef locations |
|     `group_names` |                  `Vector{String}` |   `Kora.TARGET_GROUPS` |                                Labels for the functional coral groups; must match the groups used to fit the supplied models |
|         `density` |     `Union{Int64, Vector{Int64}}` |                     20 | Maximum colony density in colonies per $\text{m}^2$; scalar applies to all locations, vector specifies per-location ceilings |
|            `area` |             `Union{Real, Vector}` |                   90.0 |                                           Reef area available for coral cover in $\text{m}^2$; scalar or per-location vector |
|          `depths` | `Union{Float64, Vector{Float64}}` |                    9.0 |                              Water depth in meters; controls bleaching mortality coefficients; scalar or per-location vector |
|   `growth_models` |           `AbstractCoralBehavior` |   `Kora.growth_models` |                                                            Fitted growth model collection, one function per functional group |
| `survival_models` |           `AbstractCoralBehavior` | `Kora.survival_models` |                                                          Fitted survival model collection, one function per functional group |


### Notes on parameter types {#Notes-on-parameter-types}

The `area` and `density` parameters each accept either a scalar or a per-location vector. When a scalar is provided, it is broadcast to all locations internally. Providing a per-location `Vector` allows you to model sites with different productive areas and different maximum standing populations.

The `group_names` parameter must match the group ordering used by the `growth_models` and `survival_models` you supply. If you use the bundled models, pass `Kora.TARGET_GROUPS` or omit the argument entirely. If you fit custom models, use the same `target_groups` argument you passed to `process_ecorrap_models`.

The `depths` parameter affects the bleaching mortality coefficient applied at each location. Shallower sites experience more intense bleaching than deeper sites for a given DHW value. The default depth of 9.0 meters is representative of the mid-shelf EcoRRAP monitoring sites.
