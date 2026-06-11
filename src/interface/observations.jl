using Missings

const BIN_ID = :surv_logclass
const TRAIN_CLASS = :class_train
const TRAIN_CLASS_MEAN_ID = :class_train_mean
const TRAIN_CLASS_STD_ID = :class_train_std
const TEST_CLASS = :class_test
const TEST_CLASS_MEAN_ID = :class_test_mean
const TEST_CLASS_STD_ID = :class_test_std

"""
    area_to_diam(area::AbstractFloat)::AbstractFloat
    area_to_diam(area::AbstractString)::Union{AbstractFloat, Missing}
    area_to_diam(_::Missing)::Missing

Convert a coral planar area (cm^2) to an equivalent circle diameter (cm).

Assumes the coral footprint is circular. The conversion solves for `d` in
`area = pi * (d/2)^2`, giving `d = sqrt(4 * area / pi)`.

String inputs are parsed to `Float64`; strings that cannot be parsed are
returned as `missing`.

# Arguments
- `area`: Planar area of the coral colony in cm^2. Accepts `AbstractFloat`,
  `AbstractString`, or `Missing`.

# Returns
- `AbstractFloat`: Equivalent circle diameter in cm when input is numeric.
- `Missing`: Returned when `area` is `missing` or is a non-numeric string.

# Examples
```jldoctest
julia> using Kora

julia> area_to_diam(Float64(pi))
2.0

julia> area_to_diam(missing)
missing
```
"""
function area_to_diam(area::AbstractFloat)::AbstractFloat
    return sqrt(4.0 * area / π)
end
function area_to_diam(area::AbstractString)::Union{AbstractFloat,Missing}
    try
        return area_to_diam(parse(Float64, area))
    catch e
        if e isa ArgumentError || e isa MethodError
            return missing
        end
        rethrow(e)
    end
end
function area_to_diam(_::Missing)::Missing
    return missing
end

"""
    adaptive_min_sample_binning(data::Vector, min_samples::Int64)::Vector{Int64}

Adaptive binning ensures each bin has at least `min_samples` samples: calculate maximum possible bins, sort data, sequentially assign samples with minimum-sample enforcement, and distribute remainders to balance bins.
"""
function adaptive_min_sample_binning(data::Vector, min_samples::Int64)::Vector{Int64}
    n = length(data)

    # Maximum bins given minimum-sample constraint:
    n_bins = n ÷ min_samples  # this is an integer division

    # Edge case: insufficient data for even one bin
    if n_bins == 0
        throw(
            ArgumentError(
                "Not enough data points for minimum bin size. Need at least $min_samples samples, got $n"
            )
        )
    end

    # n_bins is now known; apply equal-width bin-assignment:
    target_size = n ÷ n_bins  # This will be >= min_samples
    remainder = n % n_bins

    # Get indices of sorted data
    sorted_indices = sortperm(data)

    # Create bin assignments for each original data point
    bin_assignments = Vector{Int64}(undef, n)

    if remainder == 0
        bin_size = target_size
    end

    # Assign samples to bins sequentially
    current_position = 1

    for bin_idx in 1:n_bins
        if remainder != 0
            # Distribute remainder samples across first few bins
            bin_size = target_size + (bin_idx <= remainder ? 1 : 0)
        end

        # Assign samples to this bin
        this_bin = sorted_indices[current_position:(current_position + bin_size - 1)]
        bin_assignments[this_bin] .= bin_idx
        current_position = (current_position + bin_size)
    end

    return bin_assignments
end

"""
    get_growth_entries(standardized_data::DataFrame)::DataFrame

Extract and prepare rows suitable for growth model fitting from a standardized
EcoRRAP demographic dataset.

Rows are kept when `growth_use == "yes"`, `survival_use == "yes"`, and
`size != "NA"`. Rows with no recorded date between observations
(`days_t1.t2` missing) are excluded by forcing their `growth_use` to "no".
Only rows with positive linear extension (i.e., net growth) are retained.

The returned DataFrame adds the following derived columns:
- `diam`: equivalent circle diameter at observation time (cm)
- `diamnext`: equivalent circle diameter at next observation (cm)
- `logdiam`: log base-2 of `diam`
- `growth`: raw area change between observations (cm^2)
- `lin_ext`: linear extension -- diameter change between observations (cm)
- `growth_rate`: annualised linear extension (cm/yr)
- `est_1yo_growth`: estimated diameter after one year at `growth_rate`

# Arguments
- `standardized_data::DataFrame`: Output of `standardize_ecorrap_data!`. Must
  contain columns `growth_use`, `survival_use`, `size`, `sizenext`, `taxa`,
  and `days_t1.t2`.

# Returns
- `DataFrame`: Filtered and augmented subset of `standardized_data` ready for
  `organize_functional_groups`.

# See Also
`standardize_ecorrap_data!`, `get_survival_entries`, `organize_functional_groups`
"""
function get_growth_entries(standardized_data::DataFrame)::DataFrame
    # Construct masks to remove unused and missing data

    # Do not use growth data marked for use with no dates between observations!
    # Materialise the column first — Parquet2-backed StringRefVectors are read-only.
    standardized_data[!, :growth_use] = Vector{Union{Missing,String}}(standardized_data[!, :growth_use])
    growth_use_check = ismissing.(standardized_data[!, Symbol("days_t1.t2")])
    standardized_data[growth_use_check, :growth_use] .= "no"

    growth_mask = standardized_data.growth_use .== "yes"
    survived_mask = standardized_data.survival_use .== "yes"
    non_missing_size_mask = standardized_data.size .!= "NA"

    # Remove missing and unused data
    growth_data::DataFrame = standardized_data[
        growth_mask .&& survived_mask .&& non_missing_size_mask, :
    ]

    # Add diameter column and diameter Next column
    @. growth_data[!, :diam] = area_to_diam(growth_data.size)
    @. growth_data[!, :diamnext] = area_to_diam(growth_data.sizenext)

    # Add log diameter column
    @. growth_data[!, :logdiam] = log(2, growth_data.diam)

    # Add growth and linear extension entries into data frame
    @. growth_data[!, :growth] = growth_data.sizenext - growth_data.size
    @. growth_data[!, :lin_ext] = growth_data.diamnext - growth_data.diam

    days_between_obs = growth_data[!, Symbol("days_t1.t2")]
    growth_data[!, :growth_rate] .=
        passmissing(/).(growth_data.lin_ext, (days_between_obs ./ 365.25))

    @. growth_data[!, :est_1yo_growth] =
        Float64.(growth_data.diam + growth_data.growth_rate)

    # Cast taxa String15 type to string type
    growth_data[!, :taxa] .= String.(growth_data.taxa)

    no_partial_mask = growth_data.lin_ext .> 0.0
    growth_data = growth_data[no_partial_mask, :]

    return growth_data
end

"""
    get_survival_entries(standardized_data::DataFrame)::DataFrame

Extract and prepare rows suitable for survival model fitting from a standardized
EcoRRAP demographic dataset.

Rows are kept when `survival_use == "yes"`. For rows where `sizenext` is
missing or zero, the `size` column is used as a fallback so that the
`diam_mort` column (size at the mortality event) is always populated.

The returned DataFrame adds the following derived columns:
- `diam`: equivalent circle diameter at observation time (cm)
- `diam_mort`: equivalent circle diameter at next observation or mortality
  event (cm), derived from `sizenext`
- `logdiam`: natural log of `diam`

# Arguments
- `standardized_data::DataFrame`: Output of `standardize_ecorrap_data!`. Must
  contain columns `survival_use`, `size`, `sizenext`, and `taxa`.

# Returns
- `DataFrame`: Filtered and augmented subset of `standardized_data` ready for
  `organize_functional_groups`.

# See Also
`standardize_ecorrap_data!`, `get_growth_entries`, `organize_functional_groups`
"""
function get_survival_entries(standardized_data::DataFrame)::DataFrame
    # Construct masks to remove unused and missing data
    for_survival = standardized_data[:, :survival_use] .== "yes"
    for_survival[ismissing.(for_survival)] .= 0

    # If data is "missing" in the sizenext column, fill with data in `size` column
    # survival predictions use `diam_mort` (based on `sizenext`)
    # Missing entries are filled from `size`
    missing_sizenext = ismissing.(standardized_data.sizenext)
    standardized_data[missing_sizenext, :sizenext] .= standardized_data[
        missing_sizenext, :size
    ]

    zero_size = standardized_data.sizenext .== 0.0
    standardized_data[zero_size, :sizenext] .= standardized_data[zero_size, :size]

    # Only keep data marked for use for survival regressions
    survival_data::DataFrame = standardized_data[for_survival, :]

    # Insert diameter column
    survival_data[!, :diam] .= area_to_diam.(survival_data.size)
    survival_data[!, :diam_mort] .= area_to_diam.(survival_data.sizenext)

    # Add log diameter column
    survival_data[!, :logdiam] .= log.(survival_data.diam)

    # Cast taxa String15 type to string type
    survival_data[!, :taxa] .= String.(survival_data.taxa)

    return survival_data
end

function standardize_ecorrap_data!(df::DataFrame)::DataFrame
    # Make all columns lowercase
    rename!(df, Dict(n => lowercase(n) for n in names(df)))

    try
        rename!(df, :area_t1_sqcm => :size, :area_t2_sqcm => :sizenext)
    catch
    end

    try
        rename!(df, :taxon => :taxa, :survival => :surv)
    catch
    end

    # Handle differences between ecorrap data and other combined datasets.
    df.cluster = lowercase.(df.cluster)
    df[df.cluster .== "offshore_northern", :cluster] .= "offshore_north"
    df[df.cluster .== "offshore_central", :cluster] .= "offshore_central"
    df[df.cluster .== "offshore_southern", :cluster] .= "offshore_south"

    return df
end

"""
    collate_functional_groups(
        target::String,
        group_map::DataFrame,
        src_gdf::GroupedDataFrame,
        cluster_name::String;
        reef_name::Union{String,Nothing}=nothing
    )::DataFrame

Collect and concatenate all demographic records belonging to a single functional
group and reef cluster.

Looks up the species codes (`Code` column) in `group_map` whose `Cscape_group`
field contains `target`, then vertically concatenates the matching sub-groups
from `src_gdf`. A `Cscape_group` column set to `target` is appended to the
result.

Throws `ArgumentError` if no valid species codes are found for `target` in the
specified cluster (or reef, when `reef_name` is supplied).

# Arguments
- `target::String`: Functional group name to collect (matched as a substring of
  `group_map.Cscape_group`).
- `group_map::DataFrame`: Species-to-functional-group mapping table. Must have
  columns `Code` (species code string) and `Cscape_group` (functional group
  label).
- `src_gdf::GroupedDataFrame`: Demographic data grouped by `(:taxa, :cluster)`
  or `(:taxa, :cluster, :reef)` -- typically the result of applying `groupby`
  to output from `get_growth_entries` or `get_survival_entries`.
- `cluster_name::String`: Reef cluster to extract (e.g., `"offshore_central"`).
- `reef_name::Union{String,Nothing}`: Optional reef name for finer filtering.
  When `nothing` (default), all reefs within the cluster are included.

# Returns
- `DataFrame`: Concatenated records for all species mapping to `target` within
  the specified cluster, with a `Cscape_group` column added.

# See Also
`organize_functional_groups`, `get_growth_entries`, `get_survival_entries`
"""
function collate_functional_groups(
    target::String,
    group_map::DataFrame,
    src_gdf::GroupedDataFrame,
    cluster_name::String;
    reef_name::Union{String,Nothing}=nothing
)
    valid_group_codes = group_map[occursin.(target, group_map.Cscape_group), :Code]

    # Filter for only group codes that exist
    if isnothing(reef_name)
        valid_codes = filter(
            code -> haskey(src_gdf, (code, cluster_name)), valid_group_codes
        )
    else
        valid_codes = filter(
            code -> haskey(src_gdf, (code, cluster_name, reef_name)), valid_group_codes
        )
    end

    if isempty(valid_codes)
        msg = "No valid codes found for target group '$(target)' in cluster '$(cluster_name)'"
        throw(ArgumentError(msg))
    end

    if isnothing(reef_name)
        res = reduce(vcat, src_gdf[(code, cluster_name)] for code in valid_codes)
    else
        res = reduce(vcat, src_gdf[(code, cluster_name, reef_name)] for code in valid_codes)
    end

    res[!, :Cscape_group] .= target

    return res
end

"""
    organize_functional_groups(
        target_groups::Vector{String},
        group_map::DataFrame,
        gdf::GroupedDataFrame,
        cluster_name::String;
        reef::Union{String,Nothing}=nothing,
        n_bins::Int=10,
        rng::AbstractRNG=Random.default_rng()
    )::OrderedDict

Build an ordered mapping from functional group name to train/test-split
demographic data for use in model fitting.

For each name in `target_groups`, calls `collate_functional_groups` to
assemble the raw records, then applies `train_test_split!` to add size-binned
train/test class columns. The result is an `OrderedDict` keyed by group name,
preserving the order of `target_groups`.

# Arguments
- `target_groups::Vector{String}`: Ordered list of functional group names to
  process (e.g., `Kora.TARGET_GROUPS`).
- `group_map::DataFrame`: Species-to-functional-group mapping table -- see
  `collate_functional_groups` for column requirements.
- `gdf::GroupedDataFrame`: Demographic data grouped by `(:taxa, :cluster)` or
  `(:taxa, :cluster, :reef)`.
- `cluster_name::String`: Reef cluster to extract (e.g., `"offshore_central"`).
- `reef::Union{String,Nothing}`: Optional reef name. When `nothing` (default),
  all reefs in the cluster are included.
- `n_bins::Int`: Target number of size bins passed to `train_test_split!`
  (default: 10).
- `rng::AbstractRNG`: Random number generator for reproducible train/test splits.

# Returns
- `OrderedDict{String, DataFrame}`: Keys are functional group names from
  `target_groups`; values are the corresponding DataFrames with train/test
  split columns added in place.

# See Also
`collate_functional_groups`, `train_test_split!`, `process_growth_models`,
`process_survival_models`
"""
function organize_functional_groups(
    target_groups::Vector{String},
    group_map::DataFrame,
    gdf::GroupedDataFrame,
    cluster_name::String;
    reef::Union{String,Nothing}=nothing,
    n_bins=10,
    rng::AbstractRNG=Random.default_rng()
)::OrderedDict
    groupings = OrderedDict(
        fg => train_test_split!(
            collate_functional_groups(fg, group_map, gdf, cluster_name; reef_name=reef),
            n_bins;
            rng
        )
        for fg in target_groups
    )

    return groupings
end

"""
    train_test_split!(df::DataFrame, n_bins::Int; site_col::Symbol=:site_code, rng::AbstractRNG=Random.default_rng())

Split coral demographic data into training and test sets using spatially blocked splitting and size-stratified binning.

To avoid pseudoreplication and spatial autocorrelation, the split is performed at the site level (defined by `site_col`).
Approximately 60% of sites are assigned to the training set and 40% to the test set.

The function then ensures that each size bin (created via `adaptive_min_sample_binning`) contains observations
in both the training and test sets. If a bin is empty in either set, it is merged with its neighbor
that has the most observations (Bin Merging). Bootstrapped means and standard deviations of survival
rates are computed for each resulting bin-set combination.

# Arguments
- `df` : Coral demographic data containing at minimum `logdiam`, `surv`, and `site_col` columns
- `n_bins` : Target number of size bins for stratification
- `site_col` : Column name used for spatial blocking (default: `:site_code`)
- `rng` : Random number generator for reproducible splitting (default: `Random.default_rng()`)

# Modifications
The DataFrame is modified in place with the following columns added:
- `BIN_ID` : Integer identifier for the size bin
- `TRAIN_CLASS` : Bin ID if observation is in training set, 0 otherwise
- `TEST_CLASS` : Bin ID if observation is in test set, 0 otherwise
- `TRAIN_CLASS_MEAN_ID` : Bootstrapped mean survival for training observations in each bin
- `TRAIN_CLASS_STD_ID` : Bootstrapped standard deviation of survival for training observations
- `TEST_CLASS_MEAN_ID` : Bootstrapped mean survival for test observations in each bin
- `TEST_CLASS_STD_ID` : Bootstrapped standard deviation of survival for test observations

# Returns
The modified DataFrame with spatially blocked train/test assignments and bootstrapped statistics.
"""
function train_test_split!(
    df,
    n_bins;
    site_col::Symbol=:site_code,
    rng::AbstractRNG=Random.default_rng()
)
    obs_per_bin = nrow(df) ÷ n_bins
    df[!, BIN_ID] = adaptive_min_sample_binning(df.logdiam, obs_per_bin)

    # Spatial split: split by site to avoid pseudoreplication
    sites = unique(df[!, site_col])
    n_sites = length(sites)

    # Initialize class columns
    df[!, TRAIN_CLASS] .= 0
    df[!, TEST_CLASS] .= 0
    df[!, TRAIN_CLASS_MEAN_ID] .= 0.0
    df[!, TEST_CLASS_MEAN_ID] .= 0.0
    df[!, TRAIN_CLASS_STD_ID] .= 0.0
    df[!, TEST_CLASS_STD_ID] .= 0.0

    # Assign observations to train/test based on site
    local train_mask, test_mask
    if n_sites >= 2
        # Clamp to [1, n_sites-1] so both sets are always non-empty
        n_train_sites = clamp(floor(Int64, n_sites * 0.6), 1, n_sites - 1)
        train_sites = sample(rng, sites, n_train_sites; replace=false)
        test_sites = setdiff(sites, train_sites)
        train_mask = df[!, site_col] .∈ Ref(train_sites)
        test_mask = df[!, site_col] .∈ Ref(test_sites)
    else
        # Only 1 site — fall back to observation-level split
        n = nrow(df)
        n_train = max(1, floor(Int64, n * 0.6))
        train_idx = sort(sample(rng, 1:n, n_train; replace=false))
        train_mask = falses(n)
        train_mask[train_idx] .= true
        test_mask = .!train_mask
    end

    df[train_mask, TRAIN_CLASS] .= df[train_mask, BIN_ID]
    df[test_mask, TEST_CLASS] .= df[test_mask, BIN_ID]

    # Bin Merging (Option 2): Ensure every bin has data in both sets
    # Empty bins are merged with the neighbor having more data
    actual_bins = sort(unique(df[!, BIN_ID]))

    # Helper to calculate and assign bootstrap stats for a set of indices
    function assign_stats!(indices, mean_col, std_col)
        if isempty(indices)
            return nothing
        end

        obs = collect(skipmissing(df[indices, :surv]))
        if isempty(obs)
            return nothing
        end

        samp_strat = BalancedSampling(1000)
        t_mean = bootstrap(mean, obs, samp_strat).t0[1]
        t_std = bootstrap(std, obs, samp_strat).t0[1]

        df[indices, mean_col] .= t_mean
        return df[indices, std_col] .= t_std
    end

    # First pass: calculate stats for current split
    for i in actual_bins
        assign_stats!(
            findall(df[!, TRAIN_CLASS] .== i),
            TRAIN_CLASS_MEAN_ID,
            TRAIN_CLASS_STD_ID
        )
        assign_stats!(
            findall(df[!, TEST_CLASS] .== i),
            TEST_CLASS_MEAN_ID,
            TEST_CLASS_STD_ID
        )
    end

    # Consolidation pass: merge empty bins
    for i in actual_bins
        train_idx = findall(df[!, TRAIN_CLASS] .== i)
        test_idx = findall(df[!, TEST_CLASS] .== i)

        if isempty(train_idx) || isempty(test_idx)
            # Find neighbor (i-1 or i+1) with most total observations
            neighbors = filter(x -> 1 <= x <= length(actual_bins), [i - 1, i + 1])
            if isempty(neighbors)
                continue
            end

            best_neighbor = first(
                sortperm(
                    [length(findall(df[!, BIN_ID] .== actual_bins[n])) for n in neighbors];
                    rev=true
                )
            )
            neighbor_bin = actual_bins[neighbors[best_neighbor]]

            # Merge this bin into the neighbor
            df[df[!, BIN_ID] .== i, BIN_ID] .= neighbor_bin

            # Update class assignments
            df[
                (df[!, BIN_ID] .== neighbor_bin) .& train_mask,
                TRAIN_CLASS
            ] .= neighbor_bin
            df[
                (df[!, BIN_ID] .== neighbor_bin) .& test_mask,
                TEST_CLASS
            ] .= neighbor_bin

            # Recalculate stats for the merged bin
            assign_stats!(
                findall(df[!, TRAIN_CLASS] .== neighbor_bin),
                TRAIN_CLASS_MEAN_ID,
                TRAIN_CLASS_STD_ID
            )
            assign_stats!(
                findall(df[!, TEST_CLASS] .== neighbor_bin),
                TEST_CLASS_MEAN_ID,
                TEST_CLASS_STD_ID
            )
        end
    end

    return df
end

"""
    scale_feature(data::Vector{<:Real})::Vector{<:Real}

Apply z-score normalization.
"""
function scale_feature(X::Vector{<:Real})::Vector{<:Real}
    X_scaled = (X .- mean(X)) ./ std(X)

    return X_scaled
end
