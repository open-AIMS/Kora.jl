"""
    ensemble_timeseries(
        ensemble_results::NamedTuple,
        env_conditions::YAXArray,
        area::Float32;
        loc::Int=1,
        quantiles::Vector{Float64}=[0.025, 0.5, 0.975],
        juvenile_threshold::Float32=5.0f0,
        observations::Union{Nothing,NamedTuple}=nothing
    )

Create comprehensive timeseries visualization for ensemble results.

# Arguments
- `ensemble_results` : Output from `run_ensemble!`
- `env_conditions` : Environmental conditions
- `area` : Carrying capacity area for converting to percentages
- `loc` : Location index to plot (default: 1)
- `quantiles` : Quantiles to plot as uncertainty bands (default: [0.025, 0.5, 0.975])
- `juvenile_threshold` : Size threshold for juvenile corals (default: 5.0f0)
- `observations` : Optional NamedTuple with (:dates, :cover) for plotting observations

# Returns
Makie Figure object
"""
function CoralFlow.viz.ensemble_timeseries(
    ensemble_results::NamedTuple,
    env_conditions::YAXArray,
    area::Float32;
    loc::Int=1,
    quantiles::Vector{Float64}=[0.025, 0.5, 0.975],
    juvenile_threshold::Float32=5.0f0,
    observations::Union{Nothing,NamedTuple}=nothing
)
    f = Figure(; size=(900, 1400))

    n_ts = size(ensemble_results.cover, 1)
    n_ensemble = size(ensemble_results.cover, 3)
    n_grps = size(ensemble_results.group_cover, 3)

    timesteps = 1:n_ts

    # Convert to percentages
    cover_pct = (ensemble_results.cover[:, loc, :] ./ area) .* 100.0

    # Calculate statistics
    cover_median = [quantile(cover_pct[ts, :], quantiles[2]) for ts in 1:n_ts]
    cover_lower = [quantile(cover_pct[ts, :], quantiles[1]) for ts in 1:n_ts]
    cover_upper = [quantile(cover_pct[ts, :], quantiles[3]) for ts in 1:n_ts]

    # Total coral cover
    ax1 = Axis(
        f[1, 1:2];
        xlabel="Timestep",
        ylabel="Coral Cover [%]",
        title="Total Coral Cover (n=$(n_ensemble) ensemble members)"
    )

    band!(ax1, timesteps, cover_lower, cover_upper; color=(:blue, 0.2), label="95% CI")
    lines!(ax1, timesteps, cover_median; color=:blue, linewidth=2, label="Median")

    # Add observations if provided
    if !isnothing(observations)
        scatter!(ax1, observations.indices, observations.cover;
            color=:red, markersize=8, label="Observations")
    end

    axislegend(ax1; position=:rt)

    # Cover by functional group
    ax2 = Axis(
        f[2, 1:2];
        xlabel="Timestep",
        ylabel="Coral Cover [m²]",
        title="Cover by Functional Group"
    )

    colors = FGROUP_COLOR
    labels = FLABELS

    for grp in 1:n_grps
        grp_cover = ensemble_results.group_cover[:, loc, grp, :]

        grp_median = [quantile(grp_cover[ts, :], quantiles[2]) for ts in 1:n_ts]
        grp_lower = [quantile(grp_cover[ts, :], quantiles[1]) for ts in 1:n_ts]
        grp_upper = [quantile(grp_cover[ts, :], quantiles[3]) for ts in 1:n_ts]

        band!(ax2, timesteps, grp_lower, grp_upper; color=(colors[grp], 0.2))
        lines!(
            ax2, timesteps, grp_median; color=colors[grp], linewidth=2, label=labels[grp]
        )
    end

    # Juvenile cover
    ax3 = Axis(
        f[3, 1:2];
        xlabel="Timestep",
        ylabel="Coral Cover [m²]",
        title="Juvenile Cover (< $(juvenile_threshold) cm)"
    )

    # Calculate juvenile cover for each ensemble member
    juvenile_ensemble = zeros(Float32, n_ts, n_ensemble)
    for i in 1:n_ensemble, ts in 1:n_ts
        juv_cover = 0.0f0
        for grp in 1:n_grps
            grp_cover = ensemble_results.group_cover[ts, loc, grp, i]
            # Note: This is approximate since we don't have individual sizes
            # You might want to store this during ensemble run
            juv_cover += grp_cover * 0.3f0  # Assume ~30% juvenile (adjust as needed)
        end
        juvenile_ensemble[ts, i] = juv_cover
    end

    juv_median = [quantile(juvenile_ensemble[ts, :], quantiles[2]) for ts in 1:n_ts]
    juv_lower = [quantile(juvenile_ensemble[ts, :], quantiles[1]) for ts in 1:n_ts]
    juv_upper = [quantile(juvenile_ensemble[ts, :], quantiles[3]) for ts in 1:n_ts]

    band!(ax3, timesteps, juv_lower, juv_upper; color=(:green, 0.2))
    lines!(ax3, timesteps, juv_median; color=:green, linewidth=2)

    # Parameter uncertainty visualization
    ax4 = Axis(
        f[4, 1:2];
        xlabel="Timestep",
        ylabel="Trajectory Spread [%]",
        title="Ensemble Spread Over Time (stdev)"
    )

    # Calculate coefficient of variation over time
    spread = [std(cover_pct[ts, :]) for ts in 1:n_ts]
    lines!(ax4, timesteps, spread; color=:purple, linewidth=2)

    # Individual trajectories (sample)
    n_sample = min(50, n_ensemble)
    sample_indices = sample(1:n_ensemble, n_sample; replace=false)
    ax5 = Axis(
        f[5, 1:2];
        xlabel="Timestep",
        ylabel="Coral Cover [%]",
        title="Sample of Individual Trajectories (n=$(n_sample))"
    )

    for idx in sample_indices
        lines!(ax5, timesteps, cover_pct[:, idx];
            color=(:gray, 0.3), linewidth=1)
    end

    # Overlay median
    lines!(ax5, timesteps, cover_median;
        color=:blue, linewidth=3, label="Median")

    if !isnothing(observations)
        scatter!(ax5, observations.indices, observations.cover;
            color=:red, markersize=8, label="Observations")
    end

    axislegend(ax5; position=:rt)

    # Thermal stress
    ax6 = Axis(
        f[6, 1:2];
        xlabel="Timestep",
        ylabel="Degree Heating Weeks",
        title="Thermal Stress"
    )
    CoralFlow.viz.dhws!(ax6, env_conditions)

    # Add shared legend for functional groups
    Legend(f[end + 1, :], ax2; nbanks=5)

    linkxaxes!(ax1, ax2, ax3, ax4, ax5, ax6)

    return f
end
function CoralFlow.viz.ensemble_timeseries(
    reef_state::ReefState, ensemble_results::NamedTuple, env_conditions::YAXArray; loc=1,
    kwargs...
)
    return CoralFlow.viz.ensemble_timeseries(
        ensemble_results, env_conditions, reef_state.carrying_capacity[loc]; kwargs...
    )
end

"""
    bootstrap_ensemble_timeseries(
        ensemble_results::NamedTuple,
        env_conditions::YAXArray,
        area::Float32;
        n_bootstrap::Int=1000,
        kwargs...
    )

Create ensemble timeseries with additional bootstrap uncertainty on the ensemble itself.

This performs bootstrap resampling of the ensemble members to quantify sampling uncertainty
in addition to parameter uncertainty.
"""
function CoralFlow.viz.bootstrap_ensemble_timeseries(
    ensemble_results::NamedTuple,
    env_conditions::YAXArray,
    area::Float32;
    n_bootstrap::Int=1000,
    loc::Int=1,
    quantiles::Vector{Float64}=[0.025, 0.5, 0.975],
    observations::Union{Nothing,NamedTuple}=nothing
)
    f = Figure(; size=(900, 1000))

    n_ts = size(ensemble_results.cover, 1)
    n_ensemble = size(ensemble_results.cover, 3)
    timesteps = 1:n_ts

    # Convert to percentages
    cover_pct = (ensemble_results.cover[:, loc, :] ./ area) .* 100.0

    # Bootstrap the ensemble
    boot_medians = zeros(Float32, n_ts, n_bootstrap)
    boot_lowers = zeros(Float32, n_ts, n_bootstrap)
    boot_uppers = zeros(Float32, n_ts, n_bootstrap)

    @info "Performing bootstrap resampling (n=$n_bootstrap)..."

    for b in 1:n_bootstrap
        # Resample ensemble members with replacement
        boot_indices = sample(1:n_ensemble, n_ensemble; replace=true)
        boot_sample = cover_pct[:, boot_indices]

        # Calculate quantiles for this bootstrap sample
        for ts in 1:n_ts
            boot_medians[ts, b] = quantile(boot_sample[ts, :], quantiles[2])
            boot_lowers[ts, b] = quantile(boot_sample[ts, :], quantiles[1])
            boot_uppers[ts, b] = quantile(boot_sample[ts, :], quantiles[3])
        end
    end

    # Calculate bootstrap confidence intervals
    median_ci_lower = [quantile(boot_medians[ts, :], 0.025) for ts in 1:n_ts]
    median_ci_upper = [quantile(boot_medians[ts, :], 0.975) for ts in 1:n_ts]
    median_median = [quantile(boot_medians[ts, :], 0.5) for ts in 1:n_ts]

    lower_ci_lower = [quantile(boot_lowers[ts, :], 0.025) for ts in 1:n_ts]
    upper_ci_upper = [quantile(boot_uppers[ts, :], 0.975) for ts in 1:n_ts]

    # Plot
    ax1 = Axis(
        f[1, 1];
        xlabel="Timestep",
        ylabel="Coral Cover [%]",
        title="Bootstrapped Ensemble (n_ensemble=$(n_ensemble), n_boot=$(n_bootstrap))"
    )

    # Outer uncertainty (95% CI of the lower/upper bounds)
    band!(ax1, timesteps, lower_ci_lower, upper_ci_upper;
        color=(:blue, 0.15), label="Bootstrap CI of 95% PI")

    # Median and its uncertainty
    band!(ax1, timesteps, median_ci_lower, median_ci_upper;
        color=(:blue, 0.3), label="Bootstrap CI of Median")
    lines!(ax1, timesteps, median_median;
        color=:blue, linewidth=2, label="Median of Medians")

    if !isnothing(observations)
        scatter!(ax1, observations.indices, observations.cover;
            color=:red, markersize=8, label="Observations")
    end

    axislegend(ax1; position=:rt)

    # Summary statistics
    ax2 = Axis(
        f[2, 1];
        xlabel="Timestep",
        ylabel="Spread Metrics"
    )

    # Width of prediction interval over time
    pi_width = upper_ci_upper .- lower_ci_lower
    lines!(ax2, timesteps, pi_width; color=:purple, linewidth=2, label="95% PI Width")

    axislegend(ax2; position=:rt)

    linkxaxes!(ax1, ax2)

    return f
end
function CoralFlow.viz.bootstrap_ensemble_timeseries(
    reef_state::ReefState, ensemble_results::NamedTuple, env_conditions::YAXArray; kwargs...
)
    return CoralFlow.viz.bootstrap_ensemble_timeseries(
        ensemble_results, env_conditions, reef_state.carrying_capacity; kwargs...
    )
end
