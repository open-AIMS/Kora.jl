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
    growth_data[!, :diam] .= area_to_diam.(growth_data.size)
    growth_data[!, :diamnext] .= area_to_diam.(growth_data.sizenext)

    # Add log diameter column
    growth_data[!, :logdiam] .= log.(2, growth_data.diam)

    # Add growth and linear extension entries into data frame
    growth_data[!, :growth] .= growth_data.sizenext .- growth_data.size
    growth_data[!, :lin_ext] .= growth_data.diamnext .- growth_data.diam

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
    train_test_split!(df::DataFrame, n_bins::Int; rng::AbstractRNG=Random.default_rng())

Split coral demographic data into training and test sets using size-stratified binning with
bootstrapped summary statistics.

The function performs adaptive binning based on log-transformed coral diameter (`logdiam`) to
ensure adequate sample sizes per bin. Within each bin, observations are randomly assigned to
training (60%) or test (40%) sets. Bootstrapped means and standard deviations of survival
rates are computed for each bin-set combination to support subsequent model fitting and
validation.

# Arguments
- `df` : Coral demographic data containing at minimum `logdiam` and `surv` columns
- `n_bins` : Target number of size bins for stratification
- `rng` : Random number generator for reproducible train/test splitting
  (default: `Random.default_rng()`)

# Modifications
The DataFrame is modified in place with the following columns added:
- `BIN_ID` : Integer identifier for the size bin (1 to n_bins)
- `TRAIN_CLASS` : Bin ID if observation is in training set, 0 otherwise
- `TEST_CLASS` : Bin ID if observation is in test set, 0 otherwise
- `TRAIN_CLASS_MEAN_ID` : Bootstrapped mean survival for training observations in each bin
- `TRAIN_CLASS_STD_ID` : Bootstrapped standard deviation of survival for training observations
- `TEST_CLASS_MEAN_ID` : Bootstrapped mean survival for test observations in each bin
- `TEST_CLASS_STD_ID` : Bootstrapped standard deviation of survival for test observations

# Returns
The modified DataFrame with train/test assignments and bootstrapped statistics.

# Details
- Binning uses `adaptive_min_sample_binning()` to ensure minimum observations per bin
- Training/test split is 60/40 without replacement
- Bootstrap estimates use 1000 balanced samples via `BalancedSampling`
- Missing survival values are excluded from bootstrap calculations

# Example
```julia
# Prepare survival data with required columns
coral_data = get_survival_entries(raw_ecorrap_data)

# Create stratified train/test split with 10 size bins
train_test_split!(coral_data, 10; rng=Random.seed!(123))

# Access training data for a specific bin
bin_3_training = coral_data[coral_data.TRAIN_CLASS .== 3, :]
```

See also: [`adaptive_min_sample_binning`](@ref), [`get_survival_entries`](@ref),
[`fit_survival_models`](@ref)
"""
function train_test_split!(df, n_bins; rng::AbstractRNG=Random.default_rng())
    obs_per_bin = nrow(df) ÷ n_bins
    df[!, BIN_ID] = adaptive_min_sample_binning(df.logdiam, obs_per_bin)

    # If it is not training data, then the class id will be 0
    # Same for test data.
    df[!, TRAIN_CLASS] .= 0
    df[!, TEST_CLASS] .= 0

    df[!, TRAIN_CLASS_MEAN_ID] .= 0.0
    df[!, TEST_CLASS_MEAN_ID] .= 0.0

    df[!, TRAIN_CLASS_STD_ID] .= 0.0
    df[!, TEST_CLASS_STD_ID] .= 0.0

    n_bins = sort(unique(df[!, BIN_ID]))
    for i in n_bins
        class_sample = findall(df[!, BIN_ID] .== i)
        n_obs = length(class_sample)
        n_train_sample = floor(Int64, n_obs * 0.6)

        # For test/train splitting, we do *not* want to sample with replacement.
        train_sample = sample(rng, class_sample, n_train_sample; replace=false)
        test_sample = setdiff(class_sample, train_sample)

        df[train_sample, TRAIN_CLASS] .= i
        df[test_sample, TEST_CLASS] .= i

        obs = collect(skipmissing(df[train_sample, :surv]))
        samp_strat = BalancedSampling(1000)
        train_mean = bootstrap(
            mean,
            obs,
            samp_strat
        )
        train_std = bootstrap(
            std,
            obs,
            samp_strat
        )
        df[train_sample, TRAIN_CLASS_MEAN_ID] .= train_mean.t0[1]
        df[train_sample, TRAIN_CLASS_STD_ID] .= train_std.t0[1]

        obs = collect(skipmissing(df[test_sample, :surv]))
        test_mean = bootstrap(
            mean,
            obs,
            samp_strat
        )
        test_std = bootstrap(
            std,
            obs,
            samp_strat
        )
        df[test_sample, TEST_CLASS_MEAN_ID] .= test_mean.t0[1]
        df[test_sample, TEST_CLASS_STD_ID] .= test_std.t0[1]
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
