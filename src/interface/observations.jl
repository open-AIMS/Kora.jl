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
    area_to_diam(_::Missing)::Missing

Calculate the diameter of a coral given its area (in cm²) by assuming it is a circle.
If data is missing, returns `missing`.
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

Adaptive binning based on minimum samples per bin.
Determines the number of bins based on ensuring each bin has at least
`min_samples_per_bin` samples.

Approach:
1. Calculate maximum possible bins: n ÷ `min_samples_per_bin`
2. Sort the data
3. Sequentially assign samples to bins, ensuring minimum sample requirement
4. Distribute any remainder samples across bins to balance them
"""
function adaptive_min_sample_binning(data::Vector, min_samples::Int64)::Vector{Int64}
    n = length(data)

    # Calculate how many bins we can create with the minimum requirement
    n_bins = n ÷ min_samples  # this is an integer division

    # Handle edge case where we can't create any bins
    if n_bins == 0
        throw(
            ArgumentError(
                "Not enough data points for minimum bin size. Need at least $min_samples samples, got $n"
            )
        )
    end

    # Now we know n_bins, we can use similar logic to your original function
    target_size = n ÷ n_bins  # This will be >= min_samples_per_bin
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

Prep data for growth modeling based on data standardized to known format.

Removes rows not related to the calculation of growth statistics and adds diameter, log
diameter and linear extension columns to the datset of coral demographics.
This method also marks each relevant row as train/test data.

Note: This method is sufficient for the dataset with indicated date 2025-05-10.

See also:
- `standardize_ecorrap_data!()`
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

Given the standardized dataframe from containing all entries of the coral demograph data,
remove rows not related to the calculation of survival statistics and add diameter and log
diameter columns. This method also marks each relevant row as train/test data.

Survival model is trained/calibrated on data found in the `sizenext`/`diam_mort` columns.
(i.e., observed size at mortality event).

See also:
- `standardize_ecorrap_data!()`
"""
function get_survival_entries(standardized_data::DataFrame)::DataFrame
    # Construct masks to remove unused and missing data
    for_survival = standardized_data[:, :survival_use] .== "yes"
    for_survival[ismissing.(for_survival)] .= 0

    # If data is "missing" in the sizenext column, fill with data in `size` column
    # as we make survival predictions based on the data in `diam_mort` (based on `sizenext`)
    # We fill entries that are "missing" with entries from `size`
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

    # Normalise use-flag columns: parquet may emit Bool (true/false) or String ("yes"/"no").
    for col in (:growth_use, :survival_use)
        col in propertynames(df) || continue
        raw = df[!, col]
        if nonmissingtype(eltype(raw)) <: Bool
            df[!, col] = [ismissing(x) ? missing : (x ? "yes" : "no") for x in raw]
        end
    end

    return df
end

"""
    collate_functional_groups(
        target::String,
        group_map::DataFrame,
        src_gdf::GroupedDataFrame,
        cluster_name::String;
        reef_name::Union{String,Nothing}=nothing
    )

# Arguments
- `target` : Target functional group to collate data for
- `group_map` : Identified mapping between species and the functional group
- `src_gdf` : Data grouped for each taxa and cluster
- `cluster_name` : Name of target cluster as defined in `src_gdf`
- `reef_name` : Name of target reef as defined in `src_gdf`

Collate data for functional groups in the indicated cluster.
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
        n_bins=10,
        rng::AbstractRNG=Random.default_rng()
    )::OrderedDict

Organize entries for each functional group of interest into an OrderedDict.
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
    # We iterate through bins and if one is empty, merge it with the neighbor that has more data
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
