module Kora

using Printf
using OrderedCollections

using Random
using Statistics, Bootstrap
using CurveFit
using Distributions, StatsBase, StatsFuns

using CSV, DataFrames, DimensionalData

using PrecompileSignatures: @precompile_signatures
using PrecompileTools: @compile_workload

function _kora_assets_dir()
    baked = pkgdir(Kora, "assets")
    baked !== nothing && isdir(baked) && return baked
    return joinpath(Sys.BINDIR, "..", "assets")
end

include("stats.jl")
include("metrics.jl")
include("corals/corals.jl")
include("reefs/Model.jl")
include("interface/observations.jl")
include("interface/regressions.jl")
include("interface/model_io.jl")
include("interface/create_models.jl")
include("viz/viz.jl")

export
    RMSE,
    R2,
    pearson,
    spearman,
    kendall,
    ALL_METRICS

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
    recruit_cover,
    group_cover,
    group_cover_timeseries,
    juvenile_cover

export
    ReefState,
    initialize_reef,
    initialize_coral_population!,
    generate_example_environment,
    generate_environment,
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

# Methods for calibration
export
    assign_scalers!,
    mean_colony_cover_m2,
    set_population!,
    run_model,
    run_model!,
    run_ensemble!

export
    load_models,
    save_models,
    get_growth_models,
    get_survival_models,
    check_model_pair_skew

global growth_models::Union{Nothing,PolyGrowthModel} = nothing
global survival_models::Union{Nothing,PolySurvivalModel} = nothing

# Auto-generate precompilation signatures
@precompile_signatures(Kora)

@compile_workload begin
    _gm = load_models(
        joinpath(_kora_assets_dir(), "models", "offshore_north_growth_models.json")
    )
    _sm = load_models(
        joinpath(_kora_assets_dir(), "models", "offshore_north_survival_models.json")
    )
    _reef = initialize_reef(;
        n_timesteps=50,
        n_locs=20,
        density=20,
        area=90.0,
        growth_models=_gm,
        survival_models=_sm
    )
    initialize_coral_population!(_reef; rng=Xoshiro(1))
    _env = generate_example_environment(50, 20; rng=Xoshiro(42))
    run_model!(_reef, _env; rng=Xoshiro(1))

    # Cover and metrics extraction entry points
    coral_cover(_reef)
    coral_cover(_reef, 1)
    group_cover(_reef)
    group_cover(_reef, 1)
    juvenile_cover(_reef, 1)
    group_cover_timeseries(_reef)

    # DimArray reconstruction path (Base.copy)
    _reef2 = Base.copy(_reef)

    # User-facing environment wrapper
    _dhw = zeros(Float32, 50, 20)
    generate_environment(_dhw)

    # Ensemble path — 16-param (base)
    _params = zeros(Float64, 16, 1)
    _params[2:6, 1] .= 0.2  # group proportions must sum to 1.0
    run_ensemble!(_reef2, _env, _params; rng=Xoshiro(1))

    # Ensemble path — extended (scalers + recruitment), 16 + n_groups + 2 = 23 rows
    _params_ext = zeros(Float64, 23, 1)
    _params_ext[2:6, 1] .= 0.2  # group proportions
    _params_ext[17:21, 1] .= 1.0  # neutral growth scalers
    run_ensemble!(_reef2, _env, _params_ext; rng=Xoshiro(1))

    # AOT bridge path: Matrix{Float32} env + default rng (TaskLocalRNG).
    # The bridge calls run_ensemble! without an explicit rng, so Random.GLOBAL_RNG
    # (TaskLocalRNG) is used.  Without this entry the rand(TaskLocalRNG, ...) methods
    # inside run_model! are never compiled and get trimmed by juliac --trim=safe.
    _dhw_mat = generate_example_dhw(50, 20)
    _params_aot = zeros(Float64, 6, 1)
    _params_aot[1, 1] = 1.0
    _params_aot[2:6, 1] .= 0.2
    run_ensemble!(_reef2, _dhw_mat, _params_aot)
end

function __init__()
    return nothing
end

function _set_models!(gm::PolyGrowthModel, sm::PolySurvivalModel)::Nothing
    global growth_models = gm
    global survival_models = sm
    return nothing
end

end  # module Kora
