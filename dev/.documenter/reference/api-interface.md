
# Model I/O API {#Model-I/O-API}
<details class='jldocstring custom-block' open>
<summary><a id='Kora.load_models' href='#Kora.load_models'><span class="jlbinding">Kora.load_models</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
load_models(filepath::String)::Union{PolyGrowthModel, PolySurvivalModel}
```


Deserialise a model collection from a versioned JSON file previously written by `save_models`.

The file must declare a supported `format_version`, a `model_kind` of either `"growth"` or `"survival"`, and a `"type"` tag for each model entry that has been registered via `register_model_type!`. An informative error is raised if any required field is missing or if an unknown type tag is encountered.

**Arguments**
- `filepath::String` : Path to the JSON model file.
  

**Returns**

`Union{PolyGrowthModel, PolySurvivalModel}` : The deserialised model collection, ready for use as the `growth_models` or `survival_models` argument to `initialize_reef`.

**Examples**

```julia
julia> using Kora

julia> path = joinpath(pkgdir(Kora), "assets", "models",
           "offshore_north_growth_models.json");

julia> m = load_models(path);

julia> typeof(m)
Kora.PolyGrowthModel
```


**See Also**

[`save_models`](/reference/api-interface#Kora.save_models), [`register_model_type!`](/reference/api-interface#Kora.register_model_type!), [`initialize_reef`](/reference/api-reef-state#Kora.initialize_reef)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/Kora.jl/blob/319707aa5d4659465eb42066e8cff7c4b1e14e2e/src/interface/model_io.jl#L236-L271" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Kora.save_models' href='#Kora.save_models'><span class="jlbinding">Kora.save_models</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
save_models(m::PolyGrowthModel, filepath::String; region::String="")::Nothing
save_models(m::PolySurvivalModel, filepath::String; region::String="")::Nothing
```


Serialise a fitted model collection to a versioned JSON file at `filepath`.

The file records model kind (`"growth"` or `"survival"`), fit timestamp, polynomial coefficients, numeric range limits, and per-metric train/test performance. It can be reloaded without information loss using `load_models`.

**Arguments**
- `m` : Fitted model collection to serialise (`PolyGrowthModel` or `PolySurvivalModel`).
  
- `filepath::String` : Destination path including filename and `.json` extension.
  
- `region::String` : Optional label stored in the file header for reference. Does not affect the model data (default: `""`).
  

**Returns**

`Nothing`

**See Also**

[`load_models`](/reference/api-interface#Kora.load_models), [`register_model_type!`](/reference/api-interface#Kora.register_model_type!), [`fit_growth_models`](/reference/api-interface#Kora.fit_growth_models), [`fit_survival_models`](/reference/api-interface#Kora.fit_survival_models)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/Kora.jl/blob/319707aa5d4659465eb42066e8cff7c4b1e14e2e/src/interface/model_io.jl#L172-L195" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Kora.register_model_type!' href='#Kora.register_model_type!'><span class="jlbinding">Kora.register_model_type!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
register_model_type!(tag::String, deserializer::Function)::Nothing
```


Register a custom model type so that `load_models` can deserialise it from JSON.

The `tag` string must match the `"type"` field written by your serialiser into each model entry. The `deserializer` receives the raw parsed JSON `Dict` for that entry and must return a callable representing the fitted model function.

The built-in types `"PolyGrowthFunction"` and `"PolySurvivalFunction"` are registered automatically when the package is loaded via `__init__`. Registration is not thread-safe; all calls must occur before concurrent model loading begins.

**Arguments**
- `tag::String` : Unique string identifier stored in the JSON file.
  
- `deserializer::Function` : Function with signature `(entry::AbstractDict) -> Function` that reconstructs a model callable from the serialised `Dict`.
  

**Returns**

`Nothing`

**See Also**

[`load_models`](/reference/api-interface#Kora.load_models), [`save_models`](/reference/api-interface#Kora.save_models)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/Kora.jl/blob/319707aa5d4659465eb42066e8cff7c4b1e14e2e/src/interface/model_io.jl#L20-L44" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Kora.fit_growth_models' href='#Kora.fit_growth_models'><span class="jlbinding">Kora.fit_growth_models</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
fit_growth_models(
    groupings::OrderedDict{String, DataFrame}; degree=2
)::PolyGrowthModel
```


Fit polynomial growth models to grouped coral data.

**Arguments**
- `groupings` : An `OrderedDict` mapping functional-group names to `DataFrame`s, each containing `diam` (coral diameter in cm) and `growth_rate` columns, with a train/test split column produced by `get_growth_entries`.
  
- `degree` : Degree of the polynomial fitted to `log(diam)` vs `growth_rate` (default: 2).
  

**Returns**

A `PolyGrowthModel` containing one fitted function per functional group and train/test performance metrics for all [`ALL_METRICS`](/reference/api-metrics#Kora.ALL_METRICS).

**See Also**

[`fit_survival_models`](/reference/api-interface#Kora.fit_survival_models), [`save_models`](/reference/api-interface#Kora.save_models)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/Kora.jl/blob/319707aa5d4659465eb42066e8cff7c4b1e14e2e/src/interface/regressions.jl#L1-L21" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Kora.fit_survival_models' href='#Kora.fit_survival_models'><span class="jlbinding">Kora.fit_survival_models</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
fit_survival_models(
    groupings::OrderedDict{String, DataFrame}; degree=2
)
```


Fit survival models to grouped coral data using logistic regression.

**Returns**

Tuple of (models, mcfadden_r2_scores, log_likelihood_scores, brier_scores)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/Kora.jl/blob/319707aa5d4659465eb42066e8cff7c4b1e14e2e/src/interface/regressions.jl#L81-L90" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Kora.collate_functional_groups' href='#Kora.collate_functional_groups'><span class="jlbinding">Kora.collate_functional_groups</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
collate_functional_groups(
    target::String,
    group_map::DataFrame,
    src_gdf::GroupedDataFrame,
    cluster_name::String;
    reef_name::Union{String,Nothing}=nothing
)::DataFrame
```


Collect and concatenate all demographic records belonging to a single functional group and reef cluster.

Looks up the species codes (`Code` column) in `group_map` whose `Cscape_group` field contains `target`, then vertically concatenates the matching sub-groups from `src_gdf`. A `Cscape_group` column set to `target` is appended to the result.

Throws `ArgumentError` if no valid species codes are found for `target` in the specified cluster (or reef, when `reef_name` is supplied).

**Arguments**
- `target::String`: Functional group name to collect (matched as a substring of `group_map.Cscape_group`).
  
- `group_map::DataFrame`: Species-to-functional-group mapping table. Must have columns `Code` (species code string) and `Cscape_group` (functional group label).
  
- `src_gdf::GroupedDataFrame`: Demographic data grouped by `(:taxa, :cluster)` or `(:taxa, :cluster, :reef)` – typically the result of applying `groupby` to output from `get_growth_entries` or `get_survival_entries`.
  
- `cluster_name::String`: Reef cluster to extract (e.g., `"offshore_central"`).
  
- `reef_name::Union{String,Nothing}`: Optional reef name for finer filtering. When `nothing` (default), all reefs within the cluster are included.
  

**Returns**
- `DataFrame`: Concatenated records for all species mapping to `target` within the specified cluster, with a `Cscape_group` column added.
  

**See Also**

`organize_functional_groups`, `get_growth_entries`, `get_survival_entries`


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/Kora.jl/blob/319707aa5d4659465eb42066e8cff7c4b1e14e2e/src/interface/observations.jl#L279-L318" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Kora.organize_functional_groups' href='#Kora.organize_functional_groups'><span class="jlbinding">Kora.organize_functional_groups</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
organize_functional_groups(
    target_groups::Vector{String},
    group_map::DataFrame,
    gdf::GroupedDataFrame,
    cluster_name::String;
    reef::Union{String,Nothing}=nothing,
    n_bins::Int=10,
    rng::AbstractRNG=Random.default_rng()
)::OrderedDict
```


Build an ordered mapping from functional group name to train/test-split demographic data for use in model fitting.

For each name in `target_groups`, calls `collate_functional_groups` to assemble the raw records, then applies `train_test_split!` to add size-binned train/test class columns. The result is an `OrderedDict` keyed by group name, preserving the order of `target_groups`.

**Arguments**
- `target_groups::Vector{String}`: Ordered list of functional group names to process (e.g., `Kora.TARGET_GROUPS`).
  
- `group_map::DataFrame`: Species-to-functional-group mapping table – see `collate_functional_groups` for column requirements.
  
- `gdf::GroupedDataFrame`: Demographic data grouped by `(:taxa, :cluster)` or `(:taxa, :cluster, :reef)`.
  
- `cluster_name::String`: Reef cluster to extract (e.g., `"offshore_central"`).
  
- `reef::Union{String,Nothing}`: Optional reef name. When `nothing` (default), all reefs in the cluster are included.
  
- `n_bins::Int`: Target number of size bins passed to `train_test_split!` (default: 10).
  
- `rng::AbstractRNG`: Random number generator for reproducible train/test splits.
  

**Returns**
- `OrderedDict{String, DataFrame}`: Keys are functional group names from `target_groups`; values are the corresponding DataFrames with train/test split columns added in place.
  

**See Also**

`collate_functional_groups`, `train_test_split!`, `process_growth_models`, `process_survival_models`


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/Kora.jl/blob/319707aa5d4659465eb42066e8cff7c4b1e14e2e/src/interface/observations.jl#L355-L396" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Kora.get_growth_entries' href='#Kora.get_growth_entries'><span class="jlbinding">Kora.get_growth_entries</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
get_growth_entries(standardized_data::DataFrame)::DataFrame
```


Extract and prepare rows suitable for growth model fitting from a standardized EcoRRAP demographic dataset.

Rows are kept when `growth_use == "yes"`, `survival_use == "yes"`, and `size != "NA"`. Rows with no recorded date between observations (`days_t1.t2` missing) are excluded by forcing their `growth_use` to "no". Only rows with positive linear extension (i.e., net growth) are retained.

The returned DataFrame adds the following derived columns:
- `diam`: equivalent circle diameter at observation time (cm)
  
- `diamnext`: equivalent circle diameter at next observation (cm)
  
- `logdiam`: log base-2 of `diam`
  
- `growth`: raw area change between observations (cm^2)
  
- `lin_ext`: linear extension – diameter change between observations (cm)
  
- `growth_rate`: annualised linear extension (cm/yr)
  
- `est_1yo_growth`: estimated diameter after one year at `growth_rate`
  

**Arguments**
- `standardized_data::DataFrame`: Output of `standardize_ecorrap_data!`. Must contain columns `growth_use`, `survival_use`, `size`, `sizenext`, `taxa`, and `days_t1.t2`.
  

**Returns**
- `DataFrame`: Filtered and augmented subset of `standardized_data` ready for `organize_functional_groups`.
  

**See Also**

`standardize_ecorrap_data!`, `get_survival_entries`, `organize_functional_groups`


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/Kora.jl/blob/319707aa5d4659465eb42066e8cff7c4b1e14e2e/src/interface/observations.jl#L120-L151" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Kora.get_survival_entries' href='#Kora.get_survival_entries'><span class="jlbinding">Kora.get_survival_entries</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
get_survival_entries(standardized_data::DataFrame)::DataFrame
```


Extract and prepare rows suitable for survival model fitting from a standardized EcoRRAP demographic dataset.

Rows are kept when `survival_use == "yes"`. For rows where `sizenext` is missing or zero, the `size` column is used as a fallback so that the `diam_mort` column (size at the mortality event) is always populated.

The returned DataFrame adds the following derived columns:
- `diam`: equivalent circle diameter at observation time (cm)
  
- `diam_mort`: equivalent circle diameter at next observation or mortality event (cm), derived from `sizenext`
  
- `logdiam`: natural log of `diam`
  

**Arguments**
- `standardized_data::DataFrame`: Output of `standardize_ecorrap_data!`. Must contain columns `survival_use`, `size`, `sizenext`, and `taxa`.
  

**Returns**
- `DataFrame`: Filtered and augmented subset of `standardized_data` ready for `organize_functional_groups`.
  

**See Also**

`standardize_ecorrap_data!`, `get_growth_entries`, `organize_functional_groups`


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/Kora.jl/blob/319707aa5d4659465eb42066e8cff7c4b1e14e2e/src/interface/observations.jl#L197-L223" target="_blank" rel="noreferrer">source</a></Badge>

</details>

