function CoralFlow.viz.animate_population(
    reef_state::ReefState,
    loc::Int64,
    grp::Int64;
    nbins=50,
    framerate=2,
    filename="size_distribution.gif"
)
    f = Figure(size=(800, 600))
    ax = Axis(f[1, 1], xlabel="Diameter [cm]", ylabel="Density")

    pop_sample = coral_population(reef_state, 1, loc, grp)
    max_val = maximum(pop_sample)
    xlims!(ax, 0.0, max_val + 20.0)

    total_time = n_timesteps(reef_state)
    obs_points = Observable(pop_sample)

    hist!(ax, obs_points, bins=nbins, normalization=:pdf, color=(:blue, 0.3))

    grp_ids = getAxis(:group, reef_state.location_scalers)
    grp_name = grp_ids[grp]
    record(f, filename, 2:total_time; framerate=framerate) do t
        try
            obs_points[] = coral_population(reef_state, t, loc, grp)
        catch
            obs_points[] = @view [0.0f0][1:1]
        end
        ax.title = "Location $(loc) $(grp_name) - Timestep $t"

    end
end

function CoralFlow.viz.animate_population(
    reef_state::ReefState,
    env_conditions::YAXArray,
    loc::Int64,
    grp::Int64;
    nbins=50,
    framerate=2,
    filename="size_distribution.gif"
)
    f = Figure(size=(800, 800))
    ax = Axis(f[1:2, 1], xlabel="Diameter [cm]", ylabel="Density")
    ax2 = Axis(f[3, 1], xlabel="Time", ylabel="DHW")

    pop_sample = coral_population(reef_state, 1, loc, grp)
    max_val = maximum(pop_sample)
    xlims!(ax, 0.0, max_val + 20.0)

    total_time = n_timesteps(reef_state)
    obs_points = Observable{Any}(pop_sample)

    hist!(ax, obs_points, bins=nbins, normalization=:pdf, color=(:blue, 0.3))

    # DHW
    CoralFlow.viz.dhws!(ax2, env_conditions)

    # Indicate current time
    current_time = Observable(2)  # Start at timestep 2
    vlines!(ax2, current_time, color=:red, linewidth=2)

    grp_ids = getAxis(:group, reef_state.location_scalers)
    grp_name = grp_ids[grp]
    record(f, filename, 2:total_time; framerate=framerate) do t
        try
            obs_points[] = coral_population(reef_state, t, loc, grp)
            current_time[] = t  # Update the vertical line position
        catch err
            if err isa ArgumentError
                obs_points[] = @view [0.0f0][1:1]
                current_time[] = t
            else
                rethrow(err)
            end
        end

        ax.title = "Location $(loc) $(grp_name) - Timestep $t"
    end
end

function CoralFlow.viz.dhws!(ax, env_conditions)
    dhws = env_conditions[:, :, At(:dhw)].data
    series!(ax, dhws', color=:viridis)
    return nothing
end

function CoralFlow.viz.thermal_tolerance!(ax, reef_state::ReefState)
    tols = dropdims(mean(reef_state.wild_dhw_tolerances.data[:, :, :, 1], dims=2), dims=2)
    colors = [:royalblue, :lightseagreen, :mediumturquoise, :coral, :sandybrown]
    labels = ["Tabular Acropora", "Corymbose Acropora",
        "Corymbose non-Acropora", "Small massives", "Large massives"]

    timesteps = 1:size(tols, 1)
    for i in size(tols, 2):-1:1
        lines!(ax, timesteps, tols[:, i] .- tols[1, i], color=colors[i], label=labels[i])
    end

    # axislegend(ax, outside=true, tellwidth=false)

    return nothing
end

# Non-bootstrapped cover functions
function group_cover(reef_state::ReefState, ts::Int64)::Vector{Float32}
    n_grp = n_groups(reef_state)
    means = zeros(Float32, n_grp)

    for grp in 1:n_grp
        for loc in 1:n_locations(reef_state)
            pop = coral_population(reef_state, ts, loc, grp)
            means[grp] += sum(cover_cm_to_m2.(pop))
        end
        means[grp] /= n_locations(reef_state)
    end

    return means
end

function group_cover_timeseries(reef_state::ReefState)::Matrix{Float32}
    n_ts = n_timesteps(reef_state)
    n_grp = n_groups(reef_state)
    covers = zeros(Float32, n_ts, n_grp)

    for ts in 1:n_ts
        covers[ts, :] = group_cover(reef_state, ts)
    end

    return covers
end

# Non-bootstrapped juvenile cover functions
function juvenile_cover(
    reef_state::ReefState,
    ts::Int64;
    juvenile_threshold::Float32=5.0f0
)::Vector{Float32}
    n_grp = n_groups(reef_state)
    n_loc = n_locations(reef_state)
    means = zeros(Float32, n_grp)

    for grp in 1:n_grp
        total_cover = 0.0f0
        for loc in 1:n_loc
            pop = coral_population(reef_state, ts, loc, grp)
            juveniles = pop[0.0f0 .< pop .< juvenile_threshold]
            if !isempty(juveniles)
                total_cover += cover_cm_to_m2(juveniles)
            end
        end
        means[grp] = total_cover / n_loc
    end

    return means
end

function juvenile_cover_timeseries(
    reef_state::ReefState;
    juvenile_threshold::Float32=5.0f0
)::Matrix{Float32}
    n_ts = n_timesteps(reef_state)
    n_grp = n_groups(reef_state)
    means = zeros(Float32, n_ts, n_grp)

    for ts in 1:n_ts
        means[ts, :] = juvenile_cover(reef_state, ts; juvenile_threshold=juvenile_threshold)
    end

    return means
end

function _plot_cover_with_ci!(ax, means::Matrix{Float32}, lower_ci::Matrix{Float32}, upper_ci::Matrix{Float32})
    colors = [:royalblue, :lightseagreen, :mediumturquoise, :coral, :sandybrown]
    labels = ["Tabular Acropora", "Corymbose Acropora",
        "Corymbose non-Acropora", "Small massives", "Large massives"]

    timesteps = 1:size(means, 1)

    for i in size(means, 2):-1:1
        band!(ax, timesteps, lower_ci[:, i], upper_ci[:, i], color=(colors[i], 0.3))
        lines!(ax, timesteps, means[:, i], color=colors[i], label=labels[i])
    end

    axislegend(ax, outside=true, tellwidth=true)
    return nothing
end

function _plot_cover!(ax, covers::Matrix{Float32})
    colors = [:royalblue, :lightseagreen, :mediumturquoise, :coral, :sandybrown]
    labels = ["Tabular Acropora", "Corymbose Acropora",
        "Corymbose non-Acropora", "Small massives", "Large massives"]

    timesteps = 1:size(covers, 1)

    for i in size(covers, 2):-1:1
        lines!(ax, timesteps, covers[:, i], color=colors[i], label=labels[i])
    end

    return nothing
end

"""
    bootstrap_cover(
        reef_state::ReefState,
        timestep::Int64,
        n_bootstrap::Int=100;
        rng::AbstractRNG=Random.GLOBAL_RNG
    )::Matrix{Float32}

Calculate bootstrapped confidence intervals for coral cover by resampling within functional groups.
Returns (mean_cover, lower_ci, upper_ci) for each functional group, (3 ⋅ τ) where τ is the
number of functional groups.
"""
function bootstrap_cover(
    reef_state::ReefState,
    ts::Int64;
    n_bootstrap::Int=100,
    rng::AbstractRNG=Random.GLOBAL_RNG
)::Matrix{Float32}
    n_grp = n_groups(reef_state)
    n_loc = n_locations(reef_state)

    # Store results for each group
    stats = zeros(Float32, 3, n_grp)

    # Bootstrap each group separately
    for grp in 1:n_grp
        coral_diams = Float32[]

        for loc in 1:n_loc
            append!(coral_diams, cover_cm_to_m2.(coral_population(reef_state, ts, loc, grp)))
        end

        if all(coral_diams .== 0.0) || isempty(coral_diams)
            stats[:, grp] .= 0.0

            continue
        end

        bs = bootstrap(
            sum,
            coral_diams,
            BalancedSampling(n_bootstrap)
        )

        # Calculate statistics
        ci = try
            first(confint(bs, BCaConfInt(0.95)))  # stat, upper CI, lower CI
        catch err
            if !(err isa InexactError)
                rethrow(err)
            end

            # Values are constant or near constant, so return the first value.
            fill(coral_diams[1], 3)
        end

        stats[:, grp] .= ci[1], ci[2], ci[3]
    end

    return stats
end

# Main visualization functions with bootstrap option
function CoralFlow.viz.coral_cover!(ax, reef_state; n_bootstrap::Int=0)
    if n_bootstrap > 0
        n_ts = n_timesteps(reef_state)

        boot_total = zeros(Float32, n_ts)
        lower_ci = zeros(Float32, n_ts)
        upper_ci = zeros(Float32, n_ts)

        # Calculate bootstrapped statistics
        for i in 1:n_ts
            boot_total[i], lower_ci[i], upper_ci[i] = sum.(eachrow(bootstrap_cover(reef_state, i; n_bootstrap)))
        end

        timesteps = 1:n_timesteps(reef_state)
        band!(ax, timesteps, vec(lower_ci), vec(upper_ci), color=(:blue, 0.3))
        lines!(ax, timesteps, vec(boot_total), color=:blue)
    else
        covers = coral_cover(reef_state)
        series!(ax, covers', color=:dense)
    end
    return nothing
end

function CoralFlow.viz.group_cover!(ax, reef_state::ReefState; n_bootstrap::Int=0)
    if n_bootstrap > 0
        # covers = group_cover_timeseries(reef_state)

        n_ts = n_timesteps(reef_state)
        n_grps = n_groups(reef_state)

        # Calculate bootstrapped statistics
        covers = zeros(Float32, n_ts, n_grps)
        lower_ci = zeros(Float32, n_ts, n_grps)
        upper_ci = zeros(Float32, n_ts, n_grps)
        for i in 1:n_ts
            c = bootstrap_cover(reef_state, i; n_bootstrap=n_bootstrap)
            covers[i, :] .= c[1, :]
            lower_ci[i, :] .= c[2, :]
            upper_ci[i, :] .= c[3, :]
        end

        _plot_cover_with_ci!(ax, covers, lower_ci, upper_ci)
    else
        covers = group_cover_timeseries(reef_state)
        _plot_cover!(ax, covers)
    end

    return nothing
end

function CoralFlow.viz.juvenile_cover!(
    ax,
    reef_state::ReefState;
    bootstrap::Bool=false,
    juvenile_threshold::Float32=5.0f0,
    n_bootstrap::Int=100
)
    if bootstrap
        means, lower_ci, upper_ci = bootstrap_juvenile_cover_timeseries(
            reef_state, n_bootstrap; juvenile_threshold=juvenile_threshold
        )
        _plot_cover_with_ci!(ax, means, lower_ci, upper_ci)
    else
        covers = juvenile_cover_timeseries(reef_state; juvenile_threshold=juvenile_threshold)
        _plot_cover!(ax, covers)
    end
    return nothing
end

function CoralFlow.viz.population_count!(
    ax,
    reef_state::ReefState
)
    n_loc = n_locations(reef_state)
    n_steps = n_timesteps(reef_state)
    pop_counts = zeros(n_steps)

    for t in 1:n_steps
        c = 0
        for loc in 1:n_loc
            c += CoralFlow.total_population(reef_state, t, loc)
        end

        pop_counts[t] = c
    end

    lines!(ax, pop_counts)

end

function CoralFlow.viz.timeseries(
    reef_state::ReefState,
    env_conditions::YAXArray;
    n_bootstrap::Int=0,
    juvenile_threshold::Float32=5.0f0
)
    bootstrap::Bool = n_bootstrap > 0

    f = Figure(size=(900, 1400))

    # Total coral cover
    ax1 = Axis(
        f[1, 1:2],
        xlabel="Timestep",
        ylabel="Coral Cover [m²]",
        title="Total Coral Cover" * (bootstrap ? " (with 95% CI)" : "")
    )
    CoralFlow.viz.coral_cover!(ax1, reef_state; n_bootstrap=n_bootstrap)

    # Cover by functional group
    ax2 = Axis(
        f[2, 1:2],
        xlabel="Timestep",
        ylabel="Coral Cover [m²]",
        title="Cover by Functional Group" * (bootstrap ? " (with 95% CI)" : "")
    )
    CoralFlow.viz.group_cover!(ax2, reef_state; n_bootstrap=n_bootstrap)

    # Juvenile cover
    ax3 = Axis(
        f[3, 1:2],
        xlabel="Timestep",
        ylabel="Coral Cover [m²]",
        title="Juvenile Cover" * (bootstrap ? " (with 95% CI)" : "")
    )
    CoralFlow.viz.juvenile_cover!(ax3, reef_state;
        n_bootstrap=n_bootstrap,
        juvenile_threshold=juvenile_threshold
    )

    # Adaptation
    ax4 = Axis(
        f[4, 1:2],
        xlabel="Timestep",
        ylabel="Mean Tolerance [DHW]",
        title="Adaptation"
    )
    CoralFlow.viz.thermal_tolerance!(ax4, reef_state)

    # Population count
    ax5 = Axis(
        f[5, 1:2],
        xlabel="Timestep",
        ylabel="Count",
        title="Colony Count"
    )
    CoralFlow.viz.population_count!(ax5, reef_state)

    # Thermal stress
    ax6 = Axis(
        f[6, 1:2],
        xlabel="Timestep",
        ylabel="Degree Heating Weeks",
        title="Thermal Stress"
    )
    CoralFlow.viz.dhws!(ax6, env_conditions)

    Legend(f[end+1, :], ax2, nbanks=5)

    linkxaxes!(ax1, ax2, ax3, ax4, ax5, ax6)

    return f
end