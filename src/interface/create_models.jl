using Random
using Statistics, StatsBase
using CSV
using DataFrames
using Kora

function _read_parquet_file(filepath::String)::DataFrame
    error("Parquet2 must be loaded to read .parquet files: `using Parquet2`")
end

function _read_arrow_file(filepath::String)::DataFrame
    error("Arrow must be loaded to read .arrow files: `using Arrow`")
end

function _read_datafile(filepath::String)::DataFrame
    if endswith(filepath, ".parquet")
        return _read_parquet_file(filepath)
    elseif endswith(filepath, ".arrow")
        return _read_arrow_file(filepath)
    else
        return CSV.read(filepath, DataFrame; missingstring=["NA", ""])
    end
end

"""
    process_growth_models(
        ecorrap_data_file::String,
        functional_group_file::String;
        region::String = "Offshore_Central",
        degree::Int = 2,
        save_model::Bool = true,
        output_dir::String = ".",
        target_groups::Union{Vector, Nothing} = nothing,
        seed::Int = 101,
        rng::Union{AbstractRNG, Nothing} = nothing
    )

Process EcoRRAP data to create growth models for coral functional groups.

# Arguments
- `ecorrap_data_file::String`: Path to the EcoRRAP CSV data file
- `functional_group_file::String`: Path to the functional group mapping CSV file
- `region::String`: Region to process (default: "Offshore_Central")
- `degree::Int`: Polynomial degree for growth models (default: 2)
- `plot_validation::Bool=true`: Whether to plot model performance or not (default: true)
- `save_model::Bool`: Whether to serialize model to disk (default: true)
- `output_dir::String`: Directory to save model (default: Kora package assets
- `target_groups::Union{Vector, Nothing}`: Target groups to use (default: Kora.TARGET_GROUPS)
- `seed::Int`: Random seed for reproducibility (default: 101)
- `rng::Union{AbstractRNG, Nothing}`: Random number generator (default: create new with seed)

# Returns
- `NamedTuple`: Contains `growth_fits`, `growth_groupings`

# Example
```julia
results = process_growth_models(
    "EcoRRAP data for IPM_250510.csv",
    "ecorrap to cscape species.csv";
    region="Offshore_Central"
)

growth_models = results.growth_fits
```
"""
function process_growth_models(
    ecorrap_data_file::String,
    functional_group_file::String;
    region::String="Offshore_Central",
    reef::Union{String,Nothing}=nothing,
    degree::Int64=2,
    plot_validation::Bool=false,
    save_model::Bool=true,
    output_dir::String=".",
    target_groups::Union{Vector,Nothing}=Kora.TARGET_GROUPS,
    n_bins::Int64=10,
    seed::Int64=101,
    rng::Union{AbstractRNG,Nothing}=nothing
)
    # Set up random number generator if not provided
    if isnothing(rng)
        rng = Random.seed!(seed)
    end

    # Load and process EcoRRAP data
    @info "Loading EcoRRAP data from: $ecorrap_data_file"
    ecorrap_data = _read_datafile(ecorrap_data_file)

    # Standardize data
    ecorrap_data = standardize_ecorrap_data!(ecorrap_data)

    @info "Extracting growth entries..."
    ecorrap_growth = get_growth_entries(ecorrap_data)

    if !isnothing(reef)
        if reef ∉ ecorrap_growth.reef
            throw(ArgumentError("Reef $reef not in provided data file"))
        end
    end

    # Group by taxa and cluster (and reef if provided)
    growth_gdf = if isnothing(reef)
        groupby(ecorrap_growth, [:taxa, :cluster])
    else
        groupby(ecorrap_growth, [:taxa, :cluster, :reef])
    end

    # Load functional group mapping
    @info "Loading functional group mapping from: $functional_group_file"
    functional_group_map = CSV.read(functional_group_file, DataFrame; missingstring="NA")
    functional_group_map.Code .= String.(functional_group_map.Code)

    # Organize functional groups
    @info "Organizing functional groups for region: $region"
    growth_groupings = organize_functional_groups(
        target_groups, functional_group_map, growth_gdf, region; reef, n_bins
    )

    # Fit models
    @info "Fitting growth models (degree=$degree)..."
    growth_fits = Kora.fit_growth_models(growth_groupings; degree=degree)

    # Generate validation plots if requested and extension is available
    if plot_validation && (length(methods(Kora.viz.growth_performance_plots)) > 0)
        @info "Generating validation plots..."
        Kora.viz.growth_performance_plots(growth_groupings, growth_fits)
    end

    # Save model if requested
    if save_model
        fn = "$(region)"
        if !isnothing(reef)
            fn *= "_$(reef)"
        end
        _save_growth_model(growth_fits, output_dir; fn="$(fn)_growth_models", region=region)
    end

    @info "Growth model processing complete!"

    return (
        growth_fits=growth_fits,
        growth_groupings=growth_groupings
    )
end

"""
    process_survival_models(
        ecorrap_data_file::String,
        functional_group_file::String;
        region::String="Offshore_Central",
        degree::Int64=2,
        save_model::Bool=true,
        output_dir::String=".",
        plot_validation::Bool=true,
        target_groups::Union{Vector, Nothing}=nothing,
        seed::Int64=101,
        rng::Union{AbstractRNG, Nothing}=nothing
    )

Process EcoRRAP data to create survival models for coral functional groups.

# Arguments
- `ecorrap_data_file::String`: Path to the EcoRRAP CSV data file
- `functional_group_file::String`: Path to the functional group mapping CSV file
- `region::String`: Region to process (default: "Offshore_Central")
- `degree::Int`: Polynomial degree for survival models (default: 2)
- `save_model::Bool`: Whether to serialize model to disk (default: true)
- `output_dir::String`: Directory to save model (default: Kora package assets)
- `plot_validation::Bool`: Whether to generate validation plots (default: true)
- `target_groups::Union{Vector, Nothing}`: Target groups to use (default: Kora.TARGET_GROUPS)
- `seed::Int64`: Random seed for reproducibility (default: 101)
- `rng::Union{AbstractRNG, Nothing}`: Random number generator (default: create new with seed)

# Returns
- `NamedTuple`: Contains `survival_fits`, `survival_groupings`

# Example
```julia
results = process_survival_models(
    "EcoRRAP data for IPM_250510.csv",
    "ecorrap to cscape species.csv";
    region = "Offshore_Central",
    plot_validation = true
)

survival_models = results.survival_fits
```
"""
function process_survival_models(
    ecorrap_data_file::String,
    functional_group_file::String;
    region::String="Offshore_Central",
    reef::Union{String,Nothing}=nothing,
    degree::Int=2,
    save_model::Bool=true,
    output_dir::String=".",
    plot_validation::Bool=true,
    target_groups::Vector{String}=Kora.TARGET_GROUPS,
    n_bins=10,
    seed::Int=101,
    rng::Union{AbstractRNG,Nothing}=nothing
)
    # Set up random number generator
    if rng === nothing
        rng = Random.seed!(seed)
    end

    # Load and process EcoRRAP data
    @info "Loading EcoRRAP data from: $ecorrap_data_file"
    ecorrap_data = _read_datafile(ecorrap_data_file)

    # Standardize data
    ecorrap_data = standardize_ecorrap_data!(ecorrap_data)

    @info "Extracting survival entries..."
    ecorrap_survival = get_survival_entries(ecorrap_data)

    # Group by taxa and cluster
    if isnothing(reef)
        survival_gdf = groupby(ecorrap_survival, [:taxa, :cluster])
    else
        survival_gdf = groupby(ecorrap_survival, [:taxa, :cluster, :reef])
    end

    # Load functional group mapping
    @info "Loading functional group mapping from: $functional_group_file"
    functional_group_map = CSV.read(
        functional_group_file, DataFrame; missingstring="NA"
    )
    functional_group_map.Code .= String.(functional_group_map.Code)

    # Organize functional groups
    @info "Organizing functional groups for region: $region"
    surv_groupings = organize_functional_groups(
        target_groups, functional_group_map, survival_gdf, region; reef, n_bins
    )

    # Fit models
    @info "Fitting survival models (degree=$degree)..."
    surv_fits = Kora.fit_survival_models(surv_groupings; degree=degree)

    # Generate validation plots if requested and extension is available
    if plot_validation && (length(methods(Kora.viz.survival_performance_plots)) > 0)
        @info "Generating validation plots..."
        Kora.viz.survival_performance_plots(
            surv_groupings, surv_fits; target_groups=target_groups
        )
    end

    # Save model if requested
    if save_model
        fn = "$(region)"
        if !isnothing(reef)
            fn *= "_$(reef)"
        end
        _save_survival_model(
            surv_fits, output_dir; fn="$(fn)_survival_models", region=region
        )
    end

    @info "Survival model processing complete!"

    return (
        survival_fits=surv_fits,
        survival_groupings=surv_groupings
    )
end

"""
    process_ecorrap_models(
        ecorrap_data_file::String,
        functional_group_file::String;
        region::String="Offshore_Central",
        reef=nothing,
        growth_degree::Int=2,
        survival_degree::Int=2,
        save_models::Bool=true,
        output_dir::String=".",
        plot_validation::Bool=false,
        target_groups::Union{Vector, Nothing}=nothing,
        seed::Int=101,
        rng::Union{AbstractRNG, Nothing}=nothing
    )

Process EcoRRAP data to create both growth and survival models for coral functional groups.
This function combines the individual model creation functions for convenience.

# Arguments
- `ecorrap_data_file::String`: Path to the EcoRRAP CSV data file
- `functional_group_file::String`: Path to the functional group mapping CSV file
- `region::String`: Region to process (default: "Offshore_Central")
- `reef`: Reef to process
- `growth_degree::Int`: Polynomial degree for growth models (default: 2)
- `survival_degree::Int`: Polynomial degree for survival models (default: 2)
- `save_models::Bool`: Whether to serialize models to disk (default: true)
- `output_dir::String`: Directory to save models (default: Kora package assets)
- `plot_validation::Bool`: Whether to generate validation plots (default: false)
- `target_groups::Union{Vector, Nothing}`: Target groups to use (default: Kora.TARGET_GROUPS)
- `seed::Int64`: Random seed for reproducibility (default: 101)
- `rng::Union{AbstractRNG, Nothing}`: Random number generator (default: create new with seed)

# Returns
- `NamedTuple`: Contains `growth_fits`, `survival_fits`, `growth_groupings`, `survival_groupings`

# Example
```julia
results = process_ecorrap_models(
    "EcoRRAP data for IPM_250510.csv",
    "ecorrap to cscape species.csv";
    region = "Offshore_Central",
    plot_validation = true
)

# Access the fitted models
growth_models = results.growth_fits
survival_models = results.survival_fits
```
"""
function process_ecorrap_models(
    ecorrap_data_file::String,
    functional_group_file::String;
    region::String="Offshore_Central",
    reef=nothing,
    growth_degree::Int=2,
    survival_degree::Int=2,
    seed::Int=101,
    save_models::Bool=true,
    output_dir::String=".",
    plot_validation::Bool=false,
    target_groups::Vector{String}=Kora.TARGET_GROUPS,
    n_bins=15,
    rng::AbstractRNG=Random.GLOBAL_RNG
)
    # Set up shared random number generator
    Random.seed!(rng, seed)

    @info "Processing both growth and survival models..."

    # Process growth models
    growth_results = process_growth_models(
        ecorrap_data_file,
        functional_group_file;
        region=region,
        reef=reef,
        degree=growth_degree,
        seed=seed,
        save_model=save_models,
        output_dir=output_dir,
        plot_validation=plot_validation,
        target_groups=target_groups,
        n_bins=n_bins,
        rng=rng
    )

    # Process survival models
    survival_results = process_survival_models(
        ecorrap_data_file, functional_group_file;
        region=region,
        reef=reef,
        degree=survival_degree,
        seed=seed,
        save_model=save_models,
        output_dir=output_dir,
        plot_validation=plot_validation,
        target_groups=target_groups,
        n_bins=n_bins,
        rng=rng
    )

    @info "All model processing complete!"

    return (
        growth_fits=growth_results.growth_fits,
        survival_fits=survival_results.survival_fits,
        growth_groupings=growth_results.growth_groupings,
        survival_groupings=survival_results.survival_groupings
    )
end

function _make_filepath(output_dir, fn)
    # Ensure output directory exists
    mkpath(output_dir)
    out_path = joinpath(output_dir, fn * ".json")

    return out_path
end

"""
Helper function to save growth model to disk.
"""
function _save_growth_model(
    growth_fits::PolyGrowthModel, output_dir::String;
    fn::String="growth_models", region::String=""
)::Nothing
    growth_path = _make_filepath(output_dir, fn)

    @info "Saving growth models to: $growth_path"
    save_models(growth_fits, growth_path; region=region)

    return nothing
end

"""
Helper function to save survival model to disk.
"""
function _save_survival_model(
    surv_fits::PolySurvivalModel, output_dir::String;
    fn::String="survival_models", region::String=""
)::Nothing
    survival_path = _make_filepath(output_dir, fn)

    @info "Saving survival models to: $survival_path"
    save_models(surv_fits, survival_path; region=region)

    return nothing
end

"""
Helper function to save both models to disk.
"""
function _save_models(
    growth_fits::PolyGrowthModel, surv_fits::PolySurvivalModel,
    output_dir::Union{String,Nothing}; region::String=""
)::Nothing
    _save_growth_model(growth_fits, output_dir; region=region)
    _save_survival_model(surv_fits, output_dir; region=region)

    return nothing
end
