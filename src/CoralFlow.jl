module CoralFlow

using Printf
using Serialization
using OrderedCollections
using PrettyTables

using Random
using Statistics, Bootstrap
using CurveFit, GLM, MLBase
using Distributions, KernelDensity, StatsBase, StatsFuns

using CSV, DataFrames, YAXArrays
using FLoops

const ASSET_DIR = pkgdir(CoralFlow, "assets")

include("stats.jl")
include("metrics.jl")
include("corals/corals.jl")
include("reefs/ReefState.jl")
include("reefs/flow_model.jl")
include("interface/observations.jl")
include("interface/regressions.jl")
include("interface/create_models.jl")
include("viz/viz.jl")

export
    truncated_standard_normal_mean,
    truncated_normal_mean,
    truncated_normal_cdf

export
    growth_models,
    survival_models,
    GrowthModel,
    SurvivalModel

export
    cover_cm_to_m2,
    cover_cm_to_m2!,
    coral_cover,
    recruit_cover

export
    ReefState,
    initialize_reef,
    initialize_coral_population!,
    generate_example_environment,
    n_timesteps,
    n_locations,
    n_groups

export
    coral_population,
    update_wild_sample!,
    update_deployed_sample!,
    deploy_corals!

# Interface methods for known datasets
export
    area_to_diam,
    get_growth_entries,
    get_survival_entries,
    collate_functional_groups,
    organize_functional_groups

# Methods to fit growth/mortality functions
export
    fit_growth_models,
    fit_survival_models

end  # module CoralFlow
