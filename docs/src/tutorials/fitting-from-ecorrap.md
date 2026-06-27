# Fitting Models from EcoRRAP Data

## Overview

Kora ships with bundled growth and survival models fitted to offshore northern Great Barrier
Reef data from the Ecological and Reef Restoration Action Program (EcoRRAP). These models
are loaded automatically at startup and are used whenever `initialize_reef` is called without
providing custom models.

If sites are in a different region, or if models constrained by a specific reef
or time period are needed, custom models can be fitted from raw EcoRRAP demographic survey data. The fitting
process produces `PolyGrowthModel` and `PolySurvivalModel` objects. Both can be serialised to
JSON and reloaded in future sessions without refitting.

The recommended entry point for most users is `process_ecorrap_models`, which handles the
full pipeline in one call. The sections below also document the individual steps for users
who need finer control.

## Required Input Files

Two CSV files are required.

### Coral demographic survey data

The first file contains paired colony observations. Each row records the same colony at two
time points. The minimum required columns are listed in the table below.

| Column | Type | Description |
|---|---|---|
| `area_t1_sqcm` | numeric | Colony area at first observation, in $\text{cm}^2$ |
| `area_t2_sqcm` | numeric | Colony area at second observation, in $\text{cm}^2$ |
| `taxon` | string | Species or taxon code |
| `survival` | integer | 1 if colony was alive at second observation, 0 if dead |
| `growth_use` | string | "yes" if the row should be used for growth model fitting |
| `survival_use` | string | "yes" if the row should be used for survival model fitting |
| `days_t1.t2` | numeric | Number of days between the two observations |

A site or location identifier column is also expected. The `cluster` column is used during
regional grouping and must be present. A `site_code` column is used for spatial train/test
splitting and should be included if available.

The function `standardize_ecorrap_data!` normalises column names automatically before any
further processing. It renames `area_t1_sqcm` to `size`, `area_t2_sqcm` to `sizenext`,
`taxon` to `taxa`, and `survival` to `surv`. It also standardises cluster name strings
to match the expected internal values. Manual column renaming is not required before
calling any of the processing functions.

### Species-to-functional-group mapping

The second file maps each taxon code to a functional group recognised by Kora. It requires
at minimum two columns.

| Column | Type | Description |
|---|---|---|
| `Code` | string | Taxon code matching the `taxon` column in the survey data |
| `Cscape_group` | string | Functional group identifier, must be one of the five values listed below |

The five valid functional group identifiers are drawn from `Kora.TARGET_GROUPS`.

| Identifier | Common name |
|---|---|
| `acro_table` | Tabular Acropora |
| `acro_corym` | Corymbose Acropora |
| `corym_non_acro` | Branching non-Acropora |
| `small_massive` | Small Massives |
| `large_massive` | Large Massives |

Any taxon code that does not appear in the mapping file, or that maps to a group not in
`Kora.TARGET_GROUPS`, will be silently excluded from the fitted models.

## Fitting Both Models in One Call

For most workflows, `process_ecorrap_models` is the recommended entry point; it handles the full pipeline in one call. Individual steps are documented below for workflows requiring finer control.

```julia
using Kora

results = process_ecorrap_models(
    "ecorrap_demographics.csv",
    "taxon_group_mapping.csv";
    region="Offshore_North",
    output_dir="my_models"
)

growth_fits   = results.growth_fits
survival_fits = results.survival_fits
```

The `region` keyword is used to label the saved JSON files and is stored in the file header
for reference. It does not affect which data are used for fitting. The files are saved as
`<region>_growth_models.json` and `<region>_survival_models.json` in `output_dir`.

The `degree` keyword sets the polynomial degree for both models (default is 2). If different degrees are needed for growth and survival, use `process_growth_models` and
`process_survival_models` separately.

## Processing Data Step by Step

Intermediate outputs can be inspected or modified by calling the individual steps directly.

### Loading and standardising the survey data

```julia
using CSV, DataFrames, Kora

ecorrap_data = CSV.read("ecorrap_demographics.csv", DataFrame; missingstring=["NA", ""])
Kora.standardize_ecorrap_data!(ecorrap_data)
```

`standardize_ecorrap_data!` modifies the DataFrame in place and returns it. After this call,
column names follow the internal convention used by all downstream functions.

### Extracting growth and survival entries

```julia
growth_data   = Kora.get_growth_entries(ecorrap_data)
survival_data = Kora.get_survival_entries(ecorrap_data)
```

`get_growth_entries` removes rows not suitable for growth regression and adds derived columns
including colony diameter, log diameter, linear extension, and annualised growth rate.
`get_survival_entries` removes rows not suitable for survival regression and adds colony
diameter at the time of the mortality event.

### Grouping by taxa and region

```julia
using DataFrames, CSV

growth_gdf   = groupby(growth_data,   [:taxa, :cluster])
survival_gdf = groupby(survival_data, [:taxa, :cluster])
```

`taxon_group_mapping.csv` is a user-supplied lookup table that maps each taxon code (as it
appears in the raw EcoRRAP dataset) to one of Kora's functional groups. It must have four
columns:

| Column | Description |
|--------|-------------|
| `Dataset` | Source dataset name, e.g. `juv_quadrat` or `photogrammetry` |
| `Code` | Taxon identifier as it appears in the raw data |
| `Updated.name` | Full species name (may be `NA`) |
| `Cscape_group` | Functional group label used by Kora, e.g. `acro_table`, `large_massive` |

A minimal example:

```
Dataset,Code,Updated.name,Cscape_group
juv_quadrat,Acropora,,acro_table and acro_corym
juv_quadrat,Porites,,large_massive
photogrammetry,Acor,Acropora corymbose,acro_corym
photogrammetry,Atab,Acropora table,acro_table
photogrammetry,Pmas,Porites massive,large_massive
```

The `Cscape_group` values must correspond to the functional group identifiers listed in
`Kora.TARGET_GROUPS`.

```julia
group_map = CSV.read("taxon_group_mapping.csv", DataFrame; missingstring="NA")
group_map.Code .= String.(group_map.Code)

growth_groupings   = Kora.organize_functional_groups(
    Kora.TARGET_GROUPS, group_map, growth_gdf, "offshore_north"
)
survival_groupings = Kora.organize_functional_groups(
    Kora.TARGET_GROUPS, group_map, survival_gdf, "offshore_north"
)
```

The cluster name passed to `organize_functional_groups` must match the normalised cluster
label produced by `standardize_ecorrap_data!`. The three cluster labels produced after
normalisation are `"offshore_north"`, `"offshore_central"`, and `"offshore_south"`.

`organize_functional_groups` returns an `OrderedDict` mapping each functional group identifier
to its corresponding DataFrame, with a train/test split already applied.

## Fitting the Models

Once the groupings are prepared, the fitting step is straightforward.

```julia
growth_fits = Kora.fit_growth_models(growth_groupings; degree=2)
survival_fits = Kora.fit_survival_models(survival_groupings; degree=2)
```

The `degree` keyword controls the degree of the polynomial fitted to log-transformed colony
diameter. The default of 2 (quadratic) is appropriate for most EcoRRAP datasets. A higher
degree may capture more complex growth curves but is more likely to overfit if the sample
size per group is small.

Both functions return model objects that carry per-group fitted functions as well as
train and test performance metrics including RMSE, R-squared, Pearson, Spearman, and
Kendall correlations.

## Saving and Loading Model Files

### Saving to JSON

```julia
Kora.save_models(growth_fits, "my_region_growth.json"; region="My Region")
Kora.save_models(survival_fits, "my_region_survival.json"; region="My Region")
```

The JSON files record the polynomial coefficients, numeric range limits, performance metrics,
a fitted-at timestamp, and the region label. The `format_version` field is set to 1 and is
checked on load; an informative error is raised if the version is not supported by the
installed Kora release. For the full field-by-field schema and an annotated example, see
[Fitted Model Files (JSON)](../reference/input-data-reference.md#section-3-fitted-model-files-json)
in the Input Data Reference.

### Loading from JSON

```julia
gm = Kora.load_models("my_region_growth.json")
sm = Kora.load_models("my_region_survival.json")
```

`load_models` inspects the `model_kind` field in the file and returns either a
`PolyGrowthModel` or `PolySurvivalModel`. The return value is ready to pass directly to
`initialize_reef`.

### Checking that growth and survival files match

If it is unclear whether two model files were fitted from the same dataset and fitting run,
use `check_model_pair_skew` to compare their fitted-at timestamps.

```julia
Kora.check_model_pair_skew("my_region_growth.json", "my_region_survival.json")
```

If the timestamps differ by more than 24 hours (the default threshold), a warning is emitted.
A large timestamp gap suggests the files may have come from separate fitting runs on different
data subsets, which would produce inconsistently parameterised model pairs. No error is
raised; the warning is advisory.

## Using Custom Models in a Simulation

Pass the loaded models to `initialize_reef` via the `growth_models` and `survival_models`
keyword arguments.

```julia
using Random, Kora

gm = Kora.load_models("my_region_growth.json")
sm = Kora.load_models("my_region_survival.json")

reef = Kora.initialize_reef(;
    n_timesteps=50,
    n_locs=20,
    area=100.0,
    density=15,
    growth_models=gm,
    survival_models=sm
)

Kora.initialize_coral_population!(reef)

environ = Kora.generate_example_environment(50, 20)
Kora.run_model!(reef, environ; rng=Random.default_rng())
```

Any call to `initialize_reef` that does not supply explicit model arguments falls back to
the package-level defaults, which are the bundled offshore northern GBR models loaded from
`assets/models/offshore_north_growth_models.json` and
`assets/models/offshore_north_survival_models.json` at package load time.
