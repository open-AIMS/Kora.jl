module MakieExt

using Random
using Distributions

using YAXArrays
using Makie
using CoralFlow
using CoralFlow.Bootstrap


function CoralFlow.viz.animate_population(
    reef_state::ReefState,
    loc::Int64,
    grp::Int64;
    nbins=50,
    framerate=2,
    filename="size_distribution.gif"
)
    f = Figure(size=(800, 600))
    ax = Axis(f[1, 1],
        xlabel="Diameter [cm]",
        ylabel="Density"
    )
    # xlims!(ax, 0.0, 210.0)
    max_val = maximum(reef_state.pop_sample.data[:, loc, grp, :])
    xlims!(ax, 0.0, max_val+20.0)

    total_time = n_timesteps(reef_state)
    obs_points = Observable(population_sample(reef_state, 1, loc, grp))

    hist!(ax, obs_points,
          bins=nbins,
          normalization=:pdf,
          color=(:blue, 0.3))

    grp_ids = getAxis(:group, reef_state.pop_sample)
    grp_name = grp_ids[grp]
    record(f, filename, 2:total_time; framerate=framerate) do t
        obs_points[] = population_sample(reef_state, t, loc, grp)
        ax.title = "Location $(loc) $(grp_name) - Timestep $t"
        # reset_limits!(ax)
    end
end

function CoralFlow.viz.coral_cover!(ax, reef_state)::Nothing
    covers = coral_cover(reef_state)
    series!(ax, covers', color=:dense)

    return nothing
end

function CoralFlow.viz.dhws!(ax, env_conditions)
    dhws = env_conditions[:, :, At(:dhw)].data
    series!(ax, dhws', color=:viridis)

    return nothing
end

"""
    group_cover!(ax, reef_state::ReefState; n_bootstrap::Int=100)::Nothing

Stacked area chart showing cover by functional group over time with confidence intervals.
"""
function CoralFlow.viz.group_cover!(
    ax,
    reef_state::ReefState;
    n_bootstrap::Int=100,
    rng::AbstractRNG=Random.GLOBAL_RNG
)::Nothing
    # Calculate bootstrapped statistics
    means, lower_ci, upper_ci = CoralFlow.viz.bootstrap_cover_timeseries(
        reef_state, n_bootstrap; rng=rng
    )

    return CoralFlow.viz.group_cover!(ax, means, lower_ci, upper_ci)
end
function CoralFlow.viz.group_cover!(
    ax,
    means,
    lower_ci,
    upper_ci
)::Nothing
    # Define colors and labels
    colors = [:royalblue, :lightseagreen, :mediumturquoise, :coral, :sandybrown]
    labels = ["Tabular Acropora", "Corymbose Acropora",
              "Corymbose non-Acropora", "Small massives", "Large massives"]

    # Plot stacked areas with confidence intervals
    n_timesteps = size(means, 1)
    timesteps = 1:n_timesteps

    # Plot from bottom to top
    for i in size(means, 2):-1:1
        # Plot confidence interval
        band!(ax, timesteps,
              lower_ci[:, i], upper_ci[:, i],
              color=(colors[i], 0.3))

        # Plot mean line
        lines!(ax, timesteps, means[:, i],
               color=colors[i], label=labels[i])
    end

    # Place legend outside right, maintain alignment with other plots
    axislegend(ax, outside=true, tellwidth=false)  # position=:right,

    return nothing
end

"""
    total_cover!(ax, reef_state::ReefState; n_bootstrap::Int=100)::Nothing

Plot total coral cover over time with confidence intervals.
"""
function CoralFlow.viz.total_cover!(
    ax,
    reef_state::ReefState;
    n_bootstrap::Int=100,
    rng::AbstractRNG=Random.GLOBAL_RNG
)::Nothing
    means, lower_ci, upper_ci = CoralFlow.viz.bootstrap_cover_timeseries(
        reef_state, n_bootstrap; rng=rng
    )

    return CoralFlow.viz.total_cover!(ax, means, lower_ci, upper_ci)
end
function CoralFlow.viz.total_cover!(
    ax,
    means,
    lower_ci,
    upper_ci
)::Nothing
    # Sum across groups for total cover
    total_mean = sum(means, dims=2)
    total_lower = sum(lower_ci, dims=2)
    total_upper = sum(upper_ci, dims=2)

    timesteps = 1:size(means, 1)

    # Plot confidence interval
    band!(ax, timesteps, vec(total_lower), vec(total_upper),
          color=(:blue, 0.3))

    # Plot mean line
    lines!(ax, timesteps, vec(total_mean),
           color=:blue, linewidth=2)

    return nothing
end

function CoralFlow.viz.thermal_tolerance!(ax, reef_state::ReefState)::Nothing
    tols = dropdims(mean(reef_state.dhw_tolerances.data[:, :, :, 1], dims=2), dims=2)

    # Define colors and labels
    colors = [:royalblue, :lightseagreen, :mediumturquoise, :coral, :sandybrown]
    labels = ["Tabular Acropora", "Corymbose Acropora",
              "Corymbose non-Acropora", "Small massives", "Large massives"]

    timesteps = 1:size(tols, 1)

    # Plot from bottom to top
    for i in size(tols, 2):-1:1
        # Plot mean line
        lines!(ax, timesteps, tols[:, i], color=colors[i], label=labels[i])
    end

    return nothing
end

"""
    bootstrap_cover(
        reef_state::ReefState,
        timestep::Int64,
        n_bootstrap::Int=100;
        rng::AbstractRNG=Random.GLOBAL_RNG
    )::Tuple{Vector{Float32}, Vector{Float32}, Vector{Float32}}

Calculate bootstrapped confidence intervals for coral cover by resampling within functional groups.
Returns (mean_cover, lower_ci, upper_ci) for each functional group.
"""

function CoralFlow.viz.bootstrap_cover(
    reef_state::ReefState,
    ts::Int64,
    n_bootstrap::Int=100;
    rng::AbstractRNG=Random.GLOBAL_RNG
)::Tuple{Vector{Float32}, Vector{Float32}, Vector{Float32}}
    n_grp = n_groups(reef_state)
    sample_size = pop_sample_size(reef_state)
    n_loc = n_locations(reef_state)

    # Store results for each group
    means = zeros(Float32, n_grp)
    lower_ci = zeros(Float32, n_grp)
    upper_ci = zeros(Float32, n_grp)

    # Cache to collect location data for a group
    loc_data = zeros(Float32, n_loc, sample_size)

    function calc_group_mean_cover(data::AbstractVector{Float32})::Float32
        total_cover = 0.0f0
        for loc_idx in eachindex(data)
            total_cover += sum(@view(data[loc_idx, :]))
        end

        return total_cover
    end

    # Bootstrap each group separately
    for grp in 1:n_grp
        for loc in 1:n_loc
            loc_data[loc, :] = cover_cm_to_m2.(population_sample(reef_state, ts, loc, grp))
        end

        _pop = loc_data[loc_data .!= 0.0]
        if length(_pop) == 1
            means[grp] = _pop[1]
            lower_ci[grp] = _pop[1]
            upper_ci[grp] = _pop[1]
            continue
        end

        bs = bootstrap(
            calc_group_mean_cover,
            _pop,
            BalancedSampling(n_bootstrap)
        )

        # Calculate statistics
        ci = try
            first(confint(bs, BCaConfInt(0.95)))  # stat, upper CI, lower CI
        catch err
            if !(err isa InexactError)
                rethrow(err)
            end

            # Values are constant or near constant, so randomly pick 3.
            rand(_pop, 3)
        end

        means[grp] = ci[1]
        lower_ci[grp] = ci[2]
        upper_ci[grp] = ci[3]
    end

    return means, lower_ci, upper_ci
end

"""
    bootstrap_cover_timeseries(
        reef_state::ReefState,
        n_bootstrap::Int=100;
        rng::AbstractRNG=Random.GLOBAL_RNG
    )::Tuple{Matrix{Float32}, Matrix{Float32}, Matrix{Float32}}

Calculate bootstrapped confidence intervals for coral cover over time.
Returns (mean_cover, lower_ci, upper_ci) matrices of shape (timesteps, groups).
"""
function CoralFlow.viz.bootstrap_cover_timeseries(
    reef_state::ReefState,
    n_bootstrap::Int=100;
    rng::AbstractRNG=Random.GLOBAL_RNG
)::Tuple{Matrix{Float32}, Matrix{Float32}, Matrix{Float32}}
    n_ts = n_timesteps(reef_state)
    n_grp = n_groups(reef_state)

    means = zeros(Float32, n_ts, n_grp)
    lower_ci = zeros(Float32, n_ts, n_grp)
    upper_ci = zeros(Float32, n_ts, n_grp)

    for ts in 1:n_ts
        means[ts, :], lower_ci[ts, :], upper_ci[ts, :] = CoralFlow.viz.bootstrap_cover(
            reef_state, ts, n_bootstrap; rng=rng
        )
    end

    return means, lower_ci, upper_ci
end

"""
    bootstrap_juvenile_cover(
        reef_state::ReefState,
        timestep::Int64,
        n_bootstrap::Int=100;
        juvenile_threshold::Float32=5.0f0,
        rng::AbstractRNG=Random.GLOBAL_RNG
    )::Tuple{Vector{Float32}, Vector{Float32}, Vector{Float32}}

Calculate bootstrapped confidence intervals for juvenile coral cover by resampling within
functional groups. Juveniles are defined as colonies with diameter < juvenile_threshold.

Returns (mean_cover, lower_ci, upper_ci) for each functional group.
"""
function CoralFlow.viz.bootstrap_juvenile_cover(
    reef_state::ReefState,
    timestep::Int64,
    n_bootstrap::Int=100;
    juvenile_threshold::Float32=5.0f0,
    rng::AbstractRNG=Random.GLOBAL_RNG
)::Tuple{Vector{Float32}, Vector{Float32}, Vector{Float32}}
    n_grp = n_groups(reef_state)
    n_loc = n_locations(reef_state)

    # Store bootstrap results for each group
    bootstrap_means = zeros(Float32, n_bootstrap, n_grp)

    # Perform bootstrap sampling for each group
    for b in 1:n_bootstrap
        for grp in 1:n_grp
            total_cover = 0.0f0
            for loc in 1:n_loc
                # Get population sample for this group at this location
                pop_sample = population_sample(reef_state, timestep, loc, grp)

                # Filter for juveniles and sample with replacement
                juveniles = pop_sample[0.0 .< pop_sample .< juvenile_threshold]
                if !isempty(juveniles)
                    sampled_sizes = sample(rng, juveniles, length(juveniles), replace=true)
                    total_cover += cover_cm_to_m2(sampled_sizes)
                end
            end

            # Store mean cover across locations
            bootstrap_means[b, grp] = total_cover / n_loc
        end
    end

    # Calculate statistics for each group
    means = vec(mean(bootstrap_means, dims=1))
    lower_ci = vec(quantile.(eachcol(bootstrap_means), 0.025))
    upper_ci = vec(quantile.(eachcol(bootstrap_means), 0.975))

    return means, lower_ci, upper_ci
end

"""
    bootstrap_juvenile_cover_timeseries(
        reef_state::ReefState,
        n_bootstrap::Int=100;
        juvenile_threshold::Float32=5.0f0,
        rng::AbstractRNG=Random.GLOBAL_RNG
    )::Tuple{Matrix{Float32}, Matrix{Float32}, Matrix{Float32}}

Calculate bootstrapped confidence intervals for juvenile coral cover over time.
Returns (mean_cover, lower_ci, upper_ci) matrices of shape (timesteps, groups).
"""
function CoralFlow.viz.bootstrap_juvenile_cover_timeseries(
    reef_state::ReefState,
    n_bootstrap::Int=100;
    juvenile_threshold::Float32=5.0f0,
    rng::AbstractRNG=Random.GLOBAL_RNG
)::Tuple{Matrix{Float32}, Matrix{Float32}, Matrix{Float32}}
    n_ts = n_timesteps(reef_state)
    n_grp = n_groups(reef_state)

    means = zeros(Float32, n_ts, n_grp)
    lower_ci = zeros(Float32, n_ts, n_grp)
    upper_ci = zeros(Float32, n_ts, n_grp)

    for ts in 1:n_ts
        means[ts, :], lower_ci[ts, :], upper_ci[ts, :] = CoralFlow.viz.bootstrap_juvenile_cover(
            reef_state, ts, n_bootstrap; juvenile_threshold=juvenile_threshold, rng=rng
        )
    end

    return means, lower_ci, upper_ci
end

"""
    juvenile_cover!(
        ax,
        reef_state::ReefState;
        juvenile_threshold::Float32=5.0f0,
        n_bootstrap::Int=100,
        rng::AbstractRNG=Random.GLOBAL_RNG
    )::Nothing

Plot juvenile coral cover by functional group over time with confidence intervals.
"""
function CoralFlow.viz.juvenile_cover!(
    ax,
    reef_state::ReefState;
    juvenile_threshold::Float32=5.0f0,
    n_bootstrap::Int=100,
    rng::AbstractRNG=Random.GLOBAL_RNG
)::Nothing
    # Calculate bootstrapped statistics
    means, lower_ci, upper_ci = CoralFlow.viz.bootstrap_juvenile_cover_timeseries(
        reef_state, n_bootstrap; juvenile_threshold=juvenile_threshold, rng=rng
    )

    # Define colors and labels (consistent with group_cover!)
    colors = [:royalblue, :lightseagreen, :mediumturquoise, :coral, :sandybrown]
    labels = ["Tabular Acropora", "Corymbose Acropora",
              "Corymbose non-Acropora", "Small massives", "Large massives"]

    # Plot stacked areas with confidence intervals
    n_timesteps = size(means, 1)
    timesteps = 1:n_timesteps

    # Calculate cumulative sums for stacking
    cum_means = cumsum(means, dims=2)
    cum_lower = cumsum(lower_ci, dims=2)
    cum_upper = cumsum(upper_ci, dims=2)

    # Plot from bottom to top
    for i in size(means, 2):-1:1
        # Plot confidence interval
        band!(ax, timesteps,
              cum_lower[:, i], cum_upper[:, i],
              color=(colors[i], 0.3))

        # Plot mean line
        lines!(ax, timesteps, cum_means[:, i],
               color=colors[i], label=labels[i])
    end

    # Place legend outside right
    axislegend(ax, outside=true, tellwidth=false)

    return nothing
end

# Update the main timeseries function to include juvenile cover
function CoralFlow.viz.timeseries(
    reef_state::ReefState,
    env_conditions::YAXArray;
    n_bootstrap::Int=100,
    juvenile_threshold::Float32=5.0f0,
    rng::AbstractRNG=Random.GLOBAL_RNG
)

    means, lower_ci, upper_ci = CoralFlow.viz.bootstrap_cover_timeseries(
        reef_state, n_bootstrap; rng=rng
    )

    f = Figure(size=(900, 1200))

    # Total coral cover over time
    ax1 = Axis(f[1, 1:2],
        xlabel="Timestep",
        ylabel="Coral Cover [m²]",
        title="Total Coral Cover (with 95% CI)"
    )
    CoralFlow.viz.total_cover!(ax1, means, lower_ci, upper_ci)

    # Cover by functional group
    ax2 = Axis(f[2, 1:2],
        xlabel="Timestep",
        ylabel="Coral Cover [m²]",
        title="Cover by Functional Group (with 95% CI)"
    )
    CoralFlow.viz.group_cover!(ax2, means, lower_ci, upper_ci)

    # Juvenile cover by functional group
    ax3 = Axis(f[3, 1:2],
        xlabel="Timestep",
        ylabel="Coral Cover [m²]",
        title="Juvenile Cover by Functional Group (with 95% CI)"
    )
    CoralFlow.viz.juvenile_cover!(ax3, reef_state;
        juvenile_threshold=juvenile_threshold,
        n_bootstrap=n_bootstrap,
        rng=rng
    )

    ax4 = Axis(f[4, 1:2],
        xlabel="Timestep",
        ylabel="Mean Tolerance [DHW]",
        title="Adaptation"
    )
    CoralFlow.viz.thermal_tolerance!(ax4, reef_state)

    # DHW over time
    ax5 = Axis(f[5, 1:2],
        xlabel="Timestep",
        ylabel="Degree Heating Weeks",
        title="Thermal Stress"
    )
    CoralFlow.viz.dhws!(ax5, env_conditions)

    linkxaxes!(ax1, ax2, ax3, ax4, ax5)

    return f
end

# function CoralFlow.viz.survival_regression(model, grp)
#     edges = CoralFlow.bin_edges()[grp, :]
#     f, ax, sp = scatter(edges, CoralFlow.survival_rates()[grp, :])
#     lines!(model.(minimum(edges):maximum(edges)))

#     return f
# end

end