const LOGCLASS_ID = :surv_logclass
const TRAIN_CLASS = :class_train
const TRAIN_CLASS_MEAN_ID = :class_train_mean
const TRAIN_CLASS_STD_ID = :class_train_std
const TEST_CLASS = :class_test
const TEST_CLASS_MEAN_ID = :class_test_mean
const TEST_CLASS_STD_ID = :class_test_std

"""
    area_to_diam(area::AbstractFloat)::AbstractFloat

Calculate the diameter of a coral given its area (in cm²) by assuming it is a circle.
"""
function area_to_diam(area::AbstractFloat)::AbstractFloat
    return sqrt(4.0 * area / π)
end

"""
    get_growth_entries(raw_data::DataFrame; rng::AbstractRNG=Random.default_rng())::DataFrame

Prep data for growth modeling.

Removes rows not related to the calculation of growth statistics and adds diameter, log
diameter and linear extension columns to the datset of coral demographics.
This method also marks each relevant row as train/test data.

Note: This method is sufficient for the dataset with indicated date 2025-05-10.
"""
function get_growth_entries(
    raw_data::DataFrame;
    rng::AbstractRNG=Random.default_rng()
)::DataFrame
    raw_data = standardize_ecorrap_data!(raw_data)

    # Construct masks to remove unused and missing data
    growth_mask = raw_data.growth_use .== "yes"
    survived_mask = raw_data.survival_use .== "yes"
    non_missing_size_mask = raw_data.size .!= "NA"

    # Remove missing and unused data
    growth_data::DataFrame = raw_data[
        growth_mask.&&survived_mask.&&non_missing_size_mask, :
    ]

    # Add diameter column and diameter Next column
    growth_data[!, :diam] .= area_to_diam.(growth_data.size)
    growth_data[!, :diamnext] .= area_to_diam.(growth_data.sizenext)

    # Add log diameter column
    growth_data[!, :logdiam] .= log.(2, growth_data.diam)

    # Add growth and linear extension entries into data frame
    growth_data[!, :growth] .= growth_data.sizenext .- growth_data.size
    growth_data[!, :lin_ext] .= growth_data.diamnext .- growth_data.diam

    # Cast taxa String15 type to string type
    growth_data[!, :taxa] .= String.(growth_data.taxa)

    no_partial_mask = growth_data.lin_ext .> 0.0
    growth_data = growth_data[no_partial_mask, :]

    logdiam_class = quantile(growth_data[:, :logdiam], 0.1:0.1:1.0)
    growth_data.growth_logclass = map(
        x -> length(logdiam_class) - argmax(findall(x .<= logdiam_class)) + 1, growth_data.logdiam
    )

    growth_data[!, TRAIN_CLASS] .= 0
    growth_data[!, TEST_CLASS] .= 0
    for i in 1:length(logdiam_class)
        class_sample = findall(growth_data.growth_logclass .== i)
        n_obs = length(class_sample)
        n_train_sample = floor(Int64, n_obs * 0.6)
        # test_sample = 1.0 - train_sample

        train_sample = sample(rng, class_sample, n_train_sample)
        test_sample = setdiff(class_sample, train_sample)
        growth_data[train_sample, TRAIN_CLASS] .= i
        growth_data[test_sample, TEST_CLASS] .= i
    end

    return growth_data
end

"""
    get_survival_entries(raw_data::DataFrame; rng::AbstractRNG=Random.default_rng())::DataFrame

Given the csv from containing all entries of the coral demograph data, remove rows not
related to the calculation of survival statistics and add diameter and log diameter columns.
This method also marks each relevant row as train/test data.
"""
function get_survival_entries(
    raw_data::DataFrame;
    rng::AbstractRNG=Random.default_rng()
)::DataFrame
    raw_data = standardize_ecorrap_data!(raw_data)

    # Construct masks to remove unused and missing data
    for_survival = raw_data[:, :survival_use] .== "yes"
    non_missing_size_mask = raw_data.size .!= "NA"

    # Remove missing and unused data
    survival_data::DataFrame = raw_data[
        for_survival.&&non_missing_size_mask, :
    ]

    # Insert diameter column
    survival_data[!, :diam] .= area_to_diam.(survival_data.size)

    # Add log diameter column
    survival_data[!, :logdiam] .= log.(survival_data.diam)

    # Cast taxa String15 type to string type
    survival_data[!, :taxa] .= String.(survival_data.taxa)

    logdiam_class = quantile(survival_data[:, :logdiam], 0.1:0.1:1.0)

    survival_data[!, LOGCLASS_ID] = map(
        x -> length(logdiam_class) - argmax(findall(x .<= logdiam_class)) + 1, survival_data.logdiam
    )

    survival_data[!, TRAIN_CLASS] .= 0
    survival_data[!, TEST_CLASS] .= 0

    survival_data[!, TRAIN_CLASS_MEAN_ID] .= 0.0
    survival_data[!, TEST_CLASS_MEAN_ID] .= 0.0

    survival_data[!, TRAIN_CLASS_STD_ID] .= 0.0
    survival_data[!, TEST_CLASS_STD_ID] .= 0.0
    for i in 1:length(logdiam_class)
        class_sample = findall(survival_data.surv_logclass .== i)
        n_obs = length(class_sample)
        n_train_sample = floor(Int64, n_obs * 0.6)
        # test_sample = 1.0 - train_sample

        train_sample = sample(rng, class_sample, n_train_sample)
        test_sample = setdiff(class_sample, train_sample)
        survival_data[train_sample, TRAIN_CLASS] .= i
        survival_data[test_sample, TEST_CLASS] .= i

        train_mean = mean(skipmissing(survival_data[train_sample, :surv]))
        test_mean = mean(skipmissing(survival_data[test_sample, :surv]))
        survival_data[train_sample, TRAIN_CLASS_MEAN_ID] .= train_mean
        survival_data[test_sample, TEST_CLASS_MEAN_ID] .= test_mean

        train_std = std(skipmissing(survival_data[train_sample, :surv]))
        test_std = std(skipmissing(survival_data[test_sample, :surv]))
        survival_data[train_sample, TRAIN_CLASS_STD_ID] .= train_std
        survival_data[test_sample, TEST_CLASS_STD_ID] .= test_std
    end

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
    df[df.cluster.=="offshore_northern", :cluster] .= "offshore_north"
    df[df.cluster.=="offshore_central", :cluster] .= "offshore_central"
    df[df.cluster.=="offshore_southern", :cluster] .= "offshore_south"


    return df
end

"""
    collate_functional_groups(
        target::String,
        group_map::DataFrame,
        src_gdf::GroupedDataFrame,
        cluster_name::String
    )

- `target` : Target functional group to collate data for
- `group_map` : Identified mapping between species and the functional group
- `src_gdf` : Data grouped for each taxa and cluster
- `cluster_name` : Name of target cluster as defined in `src_gdf`

Collate data for functional groups in the indicated cluster.
"""
function collate_functional_groups(
    target::String,
    group_map::DataFrame,
    src_gdf::GroupedDataFrame,
    cluster_name::String
)
    valid_group_codes = group_map[occursin.(target, group_map.Cscape_group), :Code]

    # Filter for only group codes that exist
    valid_codes = filter(code -> haskey(src_gdf, (code, cluster_name)), valid_group_codes)
    if isempty(valid_codes)
        msg = "No valid codes found for target group '$(target)' in cluster '$(cluster_name)'"
        throw(ArgumentError(msg))
    end

    return reduce(vcat, src_gdf[(code, cluster_name)] for code in valid_codes)
end

"""
    organize_functional_groups(
        target_groups::Vector{String},
        group_map::DataFrame,
        gdf::GroupedDataFrame,
        cluster_name::String
    )::OrderedDict

Organize entries for each functional group of interest into an OrderedDict.
"""
function organize_functional_groups(
    target_groups::Vector{String},
    group_map::DataFrame,
    gdf::GroupedDataFrame,
    cluster_name::String
)::OrderedDict
    groupings = OrderedDict(
        fg => collate_functional_groups(fg, group_map, gdf, cluster_name)
        for fg in target_groups
    )

    return groupings
end
