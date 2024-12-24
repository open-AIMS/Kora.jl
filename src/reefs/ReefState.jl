using Distributions, KernelDensity, Random, StatsBase

using CurveFit

using YAXArrays


struct ReefState{F<:AbstractFloat,   # For Float32
                 P<:Function,        # For regression models
                 Y3<:YAXArray{F,3},  # For 3D YAXArrays
                 Y4a<:YAXArray{F,4},
                 Y4b<:YAXArray{F,4},
                 Y4c<:YAXArray{F,4}}  # For 4D YAXArrays
    pop_sample::Y4a
    growth_models::Vector{P}
    survival_models::Vector{P}
    location_scalers::Y3
    depths::Vector{F}
    dhw_tolerances::Y4b
    mortalities::Y4c
    carrying_capacity::Vector{F} # Area that supports coral in m^2
    _pop_cache::Matrix{F}        # Temporary store of population
    _growth_cache::Matrix{F}     # Temporary store of population growth
    _location_cache::Vector{F}   # Temporary store of population for a given location
    _pop_buffer::Matrix{F}       # Oversized buffer for updated population
    _location_buffer::Vector{F}  # Continually updated store of cover for all locations
    _recruit_buffer::Vector{F}   # Continually refreshed store of coral recruits
end

@inline function n_timesteps(reef_state::ReefState)::Int64
    return size(reef_state.pop_sample, :timestep)
end

@inline function n_locations(reef_state::ReefState)::Int64
    return size(reef_state.pop_sample, :location)
end

@inline function n_groups(reef_state::ReefState)::Int64
    return size(reef_state.pop_sample, :group)
end

@inline function pop_sample_size(reef_state::ReefState)::Int64
    return size(reef_state.pop_sample, :sample)
end

"""
    population_sample(reef_state::ReefState, ts::Int64, loc::Int64)::SubArray{Float32}
    population_sample(reef_state::ReefState, ts::Int64, loc::Int64, grp::Int64)::SubArray{Float32}

Retrieve sample of population for a given time and location.
"""
@inline function population_sample(reef_state::ReefState, ts::Int64, loc::Int64)::SubArray{Float32}
    return @inbounds @views reef_state.pop_sample.data[ts, loc, :, :]
end
@inline function population_sample(reef_state::ReefState, ts::Int64, loc::Int64, grp::Int64)::SubArray{Float32}
    return @inbounds @views reef_state.pop_sample.data[ts, loc, grp, :]
end

function update_sample!(reef_state::ReefState, ts::Int64, loc::Int64, group::Int64, size_dist::AbstractVector{Float32})::Nothing
    @inbounds reef_state.pop_sample.data[ts, loc, group, :] .= size_dist

    return nothing
end

function update_dhw_tol_mean!(reef_state, ts::Int64, grp::Int64, vals::AbstractVector{Float32})::Nothing
    reef_state.dhw_tolerances.data[ts, :, grp, 1] = vals

    return nothing
end
function update_dhw_tol_mean!(reef_state, ts::Int64, loc::Int64, grp::Int64, val::Float32)::Nothing
    reef_state.dhw_tolerances.data[ts, loc, grp, 1] = val

    return nothing
end

function initialize_reef(;
    n_timesteps=75, n_locs=100, n_groups=5, sample_size=1000, area=90.0, depths::Union{Float64, Vector{Float64}}=7.0,
    growth_models=CoralFlow.growth_models, survival_models=CoralFlow.survival_models
)
    group_names = Symbol.(["group$(i)" for i in 1:n_groups])

    loc_ax = (
        Dim{:scaler}([:growth, :mortality]),
        Dim{:location}(1:n_locs),
        Dim{:group}(group_names)
    )
    location_scalers = YAXArray(loc_ax, fill(1.0f0, 2, n_locs, n_groups))

    dist_ax = (
        Dim{:timestep}(1:n_timesteps),
        Dim{:location}(1:n_locs),
        Dim{:group}(group_names),
        Dim{:sample}(1:sample_size)
    )
    size_dist = YAXArray(dist_ax, zeros(Float32, n_timesteps, n_locs, n_groups, sample_size))

    dhw_ax = (
            Dim{:timestep}(1:n_timesteps),
            Dim{:location}(1:n_locs),
            Dim{:group}(group_names),
            Dim{:factor}([:mean, :stdev]),
        )
    dhw_tol = YAXArray(dhw_ax, zeros(Float32, n_timesteps, n_locs, n_groups, 2))

    mort_names = [:dhw, :cyclone]
    mort_ax = (
        Dim{:timestep}(1:n_timesteps),
        Dim{:location}(1:n_locs),
        Dim{:group}(group_names),
        Dim{:mortality}(mort_names)
    )
    mort = YAXArray(mort_ax, zeros(Float32, n_timesteps, n_locs, n_groups, length(mort_names)))

    if !(depths isa Vector)
        depths = fill(depths, n_locs)
    end

    # _pop_cache::Matrix{F}        # Temporary store of population
    # _growth_cache::Matrix{F}     # Temporary store of population growth
    # _location_cache::Vector{F}   # Temporary store of population for a given location
    # _pop_buffer::Matrix{F}       # Continually updated store of cover for all locations
    # _location_buffer::Vector{F}  # Continually updated store of cover for all locations
    # _recruit_buffer::Vector{F}   # Continually refreshed store of coral recruits
    if area isa Real
        area = fill(Float32(area), n_locs)
    elseif !(area isa Vector{<:Real})
        throw(ArgumentError("Area must be a single value or a vector indicating area of each location"))
    end
    reef_state = ReefState(
        size_dist,
        growth_models,
        survival_models,
        location_scalers,
        Float32.(depths),
        dhw_tol,
        mort,
        area,
        zeros(Float32, n_locs, sample_size),
        zeros(Float32, n_locs, sample_size),
        zeros(Float32, n_locs),
        zeros(Float32, n_locs, sample_size*2),
        zeros(Float32, n_locs),
        zeros(Float32, n_locs)
    )

    return reef_state
end

# function initialize_coral_population!(
#     reef_state::ReefState,
#     μ::Float32,
#     σ::Float32,
#     loc::Int64,
#     target_pop_size::Int64;
#     rng::AbstractRNG=Random.GLOBAL_RNG
# )::Nothing

#     # Proportions and size distributions were determined by:
#     # 1. Running the model without disturbances and selecting a time point where
#     #    trajectories appear to stabilize (t=16)
#     # 2. Assessing, for a single location, the total cover for each functional group
#     #    ```
#     #    x = sum(cover_cm_to_m2(reef_state.pop_sample[16, 1, :, :].data); dims=2)
#     #    x ./ sum(x)
#     #    ```
#     group_proportion = [
#         0.044225708,
#         0.030400556,
#         0.018606702,
#         0.4809317,
#         0.42583534,
#     ]

#     # Size distributions were then determined by fitting a LogNormal distribution to each
#     # functional group.
#     #
#     # ```
#     # fit(LogNormal, reef_state.pop_sample[16, 1, :, :].data[1, :])
#     # fit(LogNormal, reef_state.pop_sample[16, 2, :, :].data[1, :])
#     # ...
#     # ```
#     size_distribution = [
#         LogNormal{Float32}(2.2382145f0, 0.74870664f0),
#         LogNormal{Float32}(2.1610258f0, 0.64033926f0),
#         LogNormal{Float32}(1.7919482f0, 0.62202924f0),
#         LogNormal{Float32}(2.0454452f0, 1.4701339f0),
#         LogNormal{Float32}(1.9826132f0, 1.4488107f0)
#     ]

#     edges = bin_edges()
#     for grp in 1:n_groups(reef_state)
#         dist = truncated(size_distribution[grp], 0.0, maximum(edges[grp, :]))
#         sample_sizes = convert.(Float32, rand(rng, dist, target_pop_size))
#         reef_state.pop_sample.data[1, loc, grp, :] .= sample_sizes
#     end

#     reef_state.dhw_tolerances[1, :, :, At(:mean)] .= 4.0f0
#     reef_state.dhw_tolerances[:, :, :, At(:stdev)] .= 2.75f0

#     return nothing
# end
function initialize_coral_population!(
    reef_state::ReefState,
    μ::Float32,
    σ::Float32,
    loc::Int64,
    target_pop_size::Int64;
    group_proportions::Vector{Float32} = [0.1f0, 0.2f0, 0.25f0, 0.2f0, 0.25f0],
    rng::AbstractRNG=Random.GLOBAL_RNG
)::Nothing

    # Verify proportions sum to 1
    if !isapprox(sum(group_proportions), 1.0f0, atol=1e-6)
        error("Group proportions must sum to 1.0")
    end

    # Keep the same size distributions as before
    size_distribution = [
        LogNormal{Float32}(2.2382145f0, 0.74870664f0),  # tabular Acropora
        LogNormal{Float32}(2.1610258f0, 0.64033926f0),  # corymbose Acropora
        LogNormal{Float32}(1.7919482f0, 0.62202924f0),  # Pocillopora + non-Acropora corymbose
        LogNormal{Float32}(2.0454452f0, 1.4701339f0),   # Small massives and encrusting
        LogNormal{Float32}(1.9826132f0, 1.4488107f0)    # Large massives
    ]

    edges = bin_edges()
    for grp in 1:n_groups(reef_state)
        dist = truncated(size_distribution[grp], 0.0, maximum(edges[grp, :]))

        # Scale the number of samples by the desired proportion
        n_samples = round(Int, target_pop_size * group_proportions[grp])
        sample_sizes = convert.(Float32, rand(rng, dist, n_samples))

        # Pad with zeros if needed to maintain target_pop_size
        if n_samples < target_pop_size
            sample_sizes = vcat(sample_sizes, zeros(Float32, target_pop_size - n_samples))
        end

        # Clear any juveniles
        # sample_sizes[sample_sizes .< 5.0] .= 0.0

        reef_state.pop_sample.data[1, loc, grp, :] .= sample_sizes
    end

    reef_state.dhw_tolerances[1, :, 1, At(:mean)] .= 3.751612251  # tabular Acropora
    reef_state.dhw_tolerances[1, :, 2, At(:mean)] .= 4.081622683  # corymbose Acropora
    reef_state.dhw_tolerances[1, :, 3, At(:mean)] .= 4.487465256  # Pocillopora + non-Acropora corymbose
    reef_state.dhw_tolerances[1, :, 4, At(:mean)] .= 6.165751937  # Small massives and encrusting
    reef_state.dhw_tolerances[1, :, 5, At(:mean)] .= 7.153507902  # Large massives

    reef_state.dhw_tolerances[:, :, 1, At(:stdev)] .= 2.904433676  # tabular Acropora
    reef_state.dhw_tolerances[:, :, 2, At(:stdev)] .= 3.159922076  # corymbose Acropora
    reef_state.dhw_tolerances[:, :, 3, At(:stdev)] .= 3.474118416  # Pocillopora + non-Acropora corymbose
    reef_state.dhw_tolerances[:, :, 4, At(:stdev)] .= 4.773419097  # Small massives and encrusting
    reef_state.dhw_tolerances[:, :, 5, At(:stdev)] .= 5.538122776  # Large massives

    return nothing
end
function initialize_coral_population!(reef_state::ReefState; μ=30.0f0, σ=15.0f0)
    n_locs = n_locations(reef_state)
    sample_size = pop_sample_size(reef_state)
    for loc in 1:n_locs
        initialize_coral_population!(reef_state, μ, σ, loc, sample_size)
    end
end

"""
    generate_example_environment(
        n_timesteps::Int64,
        n_locations::Int64;
        rng::AbstractRNG = Random.GLOBAL_RNG,
        start_year::Int64 = 2020,
        with_dhw=true,
        warming_rate::Float32 = 0.05f0,
        seasonal_amplitude::Float32 = 1.0f0,
        dhw_threshold::Float32 = 2.0f0,
        noise_amplitude::Float32 = 0.75f0
    )
"""
function generate_example_environment(
    n_timesteps::Int64,
    n_locations::Int64;
    rng::AbstractRNG = Random.GLOBAL_RNG,
    start_year::Int64 = 2020,
    with_dhw=true,
    warming_rate::Float32 = 0.15f0,  # Increased to match observed trends
    seasonal_amplitude::Float32 = 1.2f0,  # Increased seasonal variation
    dhw_threshold::Float32 = 4.0f0,
    noise_amplitude::Float32 = 0.9f0  # Increased for more variability
)
    # Initialize arrays
    dhw_data = zeros(Float32, n_timesteps, n_locations)
    cyclone_data = zeros(Int32, n_timesteps, n_locations)

    # Calculate baseline parameters
    years = range(start_year, length=n_timesteps) .- start_year

    for loc in 1:n_locations
        # Location-specific random offset to create spatial variation
        spatial_offset = rand(rng) * 0.8f0  # Increased spatial variation

        # Track consecutive weeks above threshold
        weeks_above_threshold = 0

        # Generate temperature anomaly time series
        for t in 1:n_timesteps
            # Long-term warming trend (RCP4.5-like)
            warming_trend = warming_rate * years[t]

            # Seasonal cycle with increasing amplitude over time
            seasonal_scaling = 1.0f0 + 0.2f0 * (years[t] / years[end])
            seasonal_cycle = seasonal_amplitude * seasonal_scaling * sin(2π * (t % 12) / 12)

            # Random weather variations with increased variability over time
            weather_scaling = 1.0f0 + 0.3f0 * (years[t] / years[end])
            weather_noise = randn(rng) * noise_amplitude * weather_scaling

            # Combined temperature anomaly
            temp_anomaly = (warming_trend +
                          seasonal_cycle +
                          weather_noise +
                          spatial_offset)

            prev_dhw = t > 1 ? dhw_data[t-1, loc] : 0.0f0

            # Convert temperature anomaly to DHW with potential acute events
            if temp_anomaly > dhw_threshold
                weeks_above_threshold += 1

                # Base DHW accumulation
                dhw_accumulation = (temp_anomaly - dhw_threshold) / 4.0f0

                # Add possibility of acute temperature spikes
                if weeks_above_threshold >= 2
                    # Chance of acute event increases with warming trend
                    acute_probability = min(0.2f0 * (1.0f0 + warming_trend), 0.4f0)

                    if rand(rng) < acute_probability
                        # Generate acute spike with magnitude increasing over time
                        time_factor = years[t] / years[end]
                        base_spike = 3.0f0 + 2.0f0 * time_factor  # Spikes get larger over time
                        spike_magnitude = rand(rng) * base_spike + 2.0f0
                        dhw_accumulation += spike_magnitude
                    end
                end

                # Calculate new DHW with modified decay
                base_dhw = prev_dhw * 0.92f0  # Slower decay during heating
                raw_dhw = base_dhw + dhw_accumulation

                # Apply soft cap with fluctuations
                if raw_dhw > 20.0f0
                    excess = raw_dhw - 20.0f0
                    damping_factor = 1.0f0 / (1.0f0 + 0.3f0 * excess)

                    # Larger fluctuations in later years
                    time_factor = years[t] / years[end]
                    max_fluctuation = 4.0f0 + 2.0f0 * time_factor
                    fluctuation = (rand(rng) - 0.5f0) * min(max_fluctuation, excess * 0.6f0)

                    dhw_data[t, loc] = 20.0f0 + (excess * damping_factor) + fluctuation
                else
                    dhw_data[t, loc] = raw_dhw
                end
            else
                # Reset counter and apply faster decay when below threshold
                weeks_above_threshold = 0
                dhw_data[t, loc] = max(0.0f0, prev_dhw * 0.7f0)
            end
        end

        # Add extreme marine heatwave events
        n_extreme_events = floor(Int, n_timesteps / 12)  # More frequent events
        for _ in 1:n_extreme_events
            event_time = rand(rng, 1:n_timesteps)
            time_progress = years[event_time] / years[end]
            event_probability = time_progress * 1.8f0  # Higher probability in later years

            if rand(rng) < event_probability
                # Duration increases with time
                base_duration = 2:5
                extra_duration = rand(rng, 0:floor(Int, 3 * time_progress))
                event_duration = rand(rng, base_duration) + extra_duration

                # Base magnitude increases with time (8 DHW to 25 DHW)
                base_magnitude = 8.0f0 + 17.0f0 * time_progress
                event_magnitude = rand(rng) * 5.0f0 + base_magnitude

                # Add the extreme event with realistic onset/decline
                for t in event_time:min(event_time + event_duration, n_timesteps)
                    relative_pos = (t - event_time) / event_duration
                    if relative_pos <= 0.3f0
                        # Rapid onset
                        scaling = relative_pos / 0.3f0
                    elseif relative_pos >= 0.7f0
                        # Gradual decline
                        scaling = 1.0f0 - ((relative_pos - 0.7f0) / 0.3f0)
                    else
                        # Peak
                        scaling = 1.0f0
                    end

                    dhw_value = event_magnitude * scaling
                    dhw_data[t, loc] = max(dhw_data[t, loc], dhw_value)
                end
            end
        end
    end

    if !with_dhw
        dhw_data = zeros(Float32, n_timesteps, n_locations)
    end

    # Combine DHW and cyclone data
    env_data = cat(dhw_data, cyclone_data, dims=3)

    v_axes = (
        Dim{:timestep}(1:n_timesteps),
        Dim{:location}(1:n_locations),
        Dim{:var}([:dhw, :cyclone_category])
    )

    return YAXArray(v_axes, env_data)
end

"""
    coral_cover(reef_state::ReefState)::Matrix{Float32}
    coral_cover(reef_state::ReefState, ts::Int64)::Matrix{Float32}
    coral_cover(reef_state::ReefState, ts::Int64, loc::Int64)::Float32

Determine coral cover.
"""
function coral_cover(reef_state::ReefState)::Matrix{Float32}
    covers = zeros(Float32, size(reef_state.pop_sample)[1:2])
    for ts in axes(covers, 1)
        @inbounds covers[ts, :] = coral_cover(reef_state, ts)
    end

    return covers
end
function coral_cover(reef_state::ReefState, ts::Int64)::Vector{Float32}
    loc_covers = reef_state._location_buffer
    Threads.@threads for loc in eachindex(loc_covers)
        @inbounds loc_covers[loc] = coral_cover(reef_state, ts, loc)
    end

    return loc_covers
end
@inline function coral_cover(reef_state::ReefState, ts::Int64, loc::Int64)::Float32
    return sum(cover_cm_to_m2(@view(reef_state.pop_sample.data[ts, loc, :, :])))
end

"""
    coral_cover(ecostate::ReefState, ts::Int64, loc::Int64, grp::Int64)::Float32

Determine coral cover for a specific functional group.
"""
@inline function coral_cover(reef_state::ReefState, ts::Int64, loc::Int64, grp::Int64)::Float32
    return cover_cm_to_m2(@view(reef_state.pop_sample.data[ts, loc, grp, :]))
end

"""
    group_cover(reef_state::ReefState)::Matrix{Float32}

Retrieve coral cover by group.
"""
function group_cover(reef_state::ReefState)::Matrix{Float32}
    tmp = sum(CoralFlow.cover_cm_to_m2.(reef_state.pop_sample).data, dims=(2,4))
    return dropdims(tmp, dims=(2,4))
end

"""
    mature_coral_cover(reef_state::ReefState, ts::Int64)::Matrix{Float32}

Calculate the cover of sexually mature corals for each location and functional group.

Returns:
- Matrix{Float32} with dimensions (n_locations, n_groups) containing mature coral cover in m²
"""
function mature_coral_cover(reef_state::ReefState, ts::Int64)::Matrix{Float32}
    n_locs = n_locations(reef_state)
    n_grps = n_groups(reef_state)
    thresholds = mature_size_thresholds()

    mature_cover = zeros(Float32, n_locs, n_grps)
    for grp in 1:n_grps
        for loc in 1:n_locs
            # Get population sample for this location/group
            pop = population_sample(reef_state, ts, loc, grp)
            # Sum cover of colonies above maturity threshold
            mature_cover[loc, grp] = cover_cm_to_m2(pop[pop .>= thresholds[grp]])
        end
    end

    return mature_cover
end

"""
Temporary mock function to determine cover of new recruits
"""
function recruit_cover(recruits::Array{Float32})::Vector{Float32}
    loc_cover = zeros(Float32, axes(recruits, 1))  # zeros for each location
    tmp_cover = cover_cm_to_m2.(recruits)
    @floop for i in eachindex(loc_cover)
        @inbounds loc_cover[i] = sum(@view(tmp_cover[i, :, :]))
    end

    return loc_cover
end

function recruit_cover(ecostate::ReefState, recruits::Array{Float32})
    loc_cover = ecostate._recruit_buffer

    tmp_cover = cover_cm_to_m2.(recruits)
    @floop for i in eachindex(loc_cover)
        @inbounds loc_cover[i] = sum(@view(tmp_cover[i, :, :]))
    end

    return loc_cover
end
