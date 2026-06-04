module Kora

using Printf
using OrderedCollections

using Random
using Statistics, Bootstrap
using CurveFit, GLM, MLBase, Interpolations
using Distributions, KernelDensity, StatsBase, StatsFuns

using CSV, DataFrames, YAXArrays

using PrecompileSignatures: @precompile_signatures
using PrecompileTools: @compile_workload

const ASSET_DIR = pkgdir(Kora, "assets")

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
    set_population!,
    run_model,
    run_model!,
    run_ensemble!

export
    load_models,
    save_models,
    get_growth_models,
    get_survival_models,
    check_model_pair_skew,
    register_model_type!

# Auto-generate precompilation signatures
@precompile_signatures(Kora)

@compile_workload begin
    register_model_type!("PolyGrowthFunction", _deserialize_poly_growth)
    register_model_type!("PolySurvivalFunction", _deserialize_poly_survival)
    _gm = load_models(joinpath(ASSET_DIR, "models", "offshore_north_growth_models.json"))
    _sm = load_models(joinpath(ASSET_DIR, "models", "offshore_north_survival_models.json"))
    _reef = initialize_reef(; n_timesteps=50, n_locs=20, density=20, area=90.0, growth_models=_gm, survival_models=_sm)
    initialize_coral_population!(_reef; rng=Xoshiro(1))
    _env = generate_example_environment(50, 20; rng=Xoshiro(42))
    run_model!(_reef, _env; rng=Xoshiro(1))
end

function __init__()
    # Populate registry at load-time, not precompile-time.
    # Functions are defined at include-time; only the dict insertion happens here.
    register_model_type!("PolyGrowthFunction", _deserialize_poly_growth)
    register_model_type!("PolySurvivalFunction", _deserialize_poly_survival)

    _growth_path = joinpath(ASSET_DIR, "models", "offshore_north_growth_models.json")
    _survival_path = joinpath(ASSET_DIR, "models", "offshore_north_survival_models.json")

    global growth_models = try
        load_models(_growth_path)
    catch e
        @warn "Pre-defined growth models could not be loaded." exception = e
        nothing
    end

    global survival_models = try
        load_models(_survival_path)
    catch e
        @warn "Pre-defined survival models could not be loaded." exception = e
        nothing
    end

    if !isnothing(growth_models) && !isnothing(survival_models)
        check_model_pair_skew(_growth_path, _survival_path)
    end
end

end  # module Kora
