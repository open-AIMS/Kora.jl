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
function Kora.viz.ensemble_timeseries(
    reef_state::ReefState,
    ensemble_results::NamedTuple,
    env_conditions::YAXArray;
    loc::Int=1,
    quantiles::Vector{Float64}=[0.025, 0.5, 0.975],
    observations::Union{Nothing,NamedTuple}=nothing
)
    area = reef_state.carrying_capacity

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
        title="Juvenile Cover"
    )

    # Calculate juvenile cover for each ensemble member
    ens_juveniles = ensemble_results.juvenile_cover
    juvenile_ensemble = zeros(Float32, n_ts, n_grps, n_ensemble)
    for i in 1:n_ensemble, ts in 1:n_ts
        for grp in 1:n_grps
            juvenile_ensemble[ts, grp, i] += ens_juveniles[ts, loc, grp, i]
        end
    end

    for grp in 1:n_grps
        grp_cover = juvenile_ensemble[:, grp, :]

        grp_median = [quantile(grp_cover[ts, :], quantiles[2]) for ts in 1:n_ts]
        grp_lower = [quantile(grp_cover[ts, :], quantiles[1]) for ts in 1:n_ts]
        grp_upper = [quantile(grp_cover[ts, :], quantiles[3]) for ts in 1:n_ts]

        band!(ax3, timesteps, grp_lower, grp_upper; color=(colors[grp], 0.2))
        lines!(
            ax3, timesteps, grp_median; color=colors[grp], linewidth=2, label=labels[grp]
        )
    end

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

    # # Individual trajectories (sample)
    # n_sample = min(50, n_ensemble)
    # sample_indices = sample(1:n_ensemble, n_sample; replace=false)
    # ax5 = Axis(
    #     f[5, 1:2];
    #     xlabel="Timestep",
    #     ylabel="Coral Cover [%]",
    #     title="Sample of Individual Trajectories (n=$(n_sample))"
    # )

    # for idx in sample_indices
    #     lines!(ax5, timesteps, cover_pct[:, idx];
    #         color=(:gray, 0.3), linewidth=1)
    # end

    # # Overlay median
    # lines!(ax5, timesteps, cover_median;
    #     color=:blue, linewidth=3, label="Median")

    # if !isnothing(observations)
    #     scatter!(ax5, observations.indices, observations.cover;
    #         color=:red, markersize=8, label="Observations")
    # end

    # axislegend(ax5; position=:rt)

    # Thermal adaptation
    ax5 = Axis(
        f[5, 1:2];
        xlabel="Timestep",
        ylabel="Mean Tolerance [DHW]",
        title="Adaptation"
    )

    # Extract tolerance data from ensemble results
    # wild_dhw_tolerances: [timestep, location, group, mean_stdev, ensemble]
    mean_tols_ensemble = ensemble_results.wild_dhw_tolerances[:, loc, :, 1, :]  # Get mean tolerance

    for grp in 1:n_grps
        # Get tolerance change relative to initial value for each ensemble member
        grp_tols = mean_tols_ensemble[:, grp, :]  # [timestep, ensemble]
        initial_tols = grp_tols[1, :]  # Initial values for each ensemble member

        # Calculate change from initial
        tol_change = zeros(Float32, n_ts, n_ensemble)
        for i in 1:n_ensemble
            tol_change[:, i] = grp_tols[:, i] .- initial_tols[i]
        end

        # Calculate statistics
        tol_median = [quantile(tol_change[ts, :], quantiles[2]) for ts in 1:n_ts]
        tol_lower = [quantile(tol_change[ts, :], quantiles[1]) for ts in 1:n_ts]
        tol_upper = [quantile(tol_change[ts, :], quantiles[3]) for ts in 1:n_ts]

        band!(
            ax5, timesteps, tol_lower, tol_upper; color=(colors[grp], 0.2)
        )
        lines!(
            ax5, timesteps, tol_median; color=colors[grp], linewidth=2,
            label=labels[grp]
        )
    end

    # Thermal stress
    ax6 = Axis(
        f[6, 1:2];
        xlabel="Timestep",
        ylabel="Degree Heating Weeks",
        title="Thermal Stress"
    )
    Kora.viz.dhws!(ax6, env_conditions)

    # Add shared legend for functional groups
    Legend(f[end + 1, :], ax2; nbanks=5)

    linkxaxes!(ax1, ax2, ax3, ax4, ax5, ax6)

    return f
end
