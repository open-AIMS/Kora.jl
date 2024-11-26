using CurveFit

# Private consts - only used in this file
const _euler_f32b = Float32(ℯ)

"""
    space_constraint(x, k; x0=0.96)

Modify growth as coral cover approaches total habitable area.

# Arguments
- `x` : proportion of available space currently taken
- `k` : steepness parameter, where the larger the value, the faster the drop (suggest k=20).
- `x0` : Inflection point where growth begins to be constrained (default: 96% of available area).
"""
@inline function space_constraint(x::F, k::F; x0::F=0.96f0)::F where {F<:Float32}
    return 1.0f0 / (1.0f0 + _euler_f32b^(k*(x - x0)))
end


# Growth models need to account for available space... (see space constraint)
@inline function growth(model::Polynomial{F}, ext::F)::F where {F<:Float32}
    if ext < 5.0f0
        return 1.0f0
    end
    if ext > 150.0f0
        return 0.1f0
    end

    return max(model(ext), 0.1f0)
end
function growth!(model::Polynomial{F}, ext::AbstractMatrix{F})::Nothing where {F<:Float32}
    Threads.@threads for i in eachindex(ext)
        # Combine model evaluation, max, and scaling in one pass
        @inbounds ext[i] += growth(model, ext[i])
    end

    return nothing
end


function update_size_distribution!(
    ecostate::ReefState,
    ts::Int64,
    loc::Int64,
    grp::Int64,
    pop::AbstractArray{Float32},
    class_diams::Matrix{Float32};
    rng::AbstractRNG=Random.GLOBAL_RNG
)::Nothing
    _pop = pop[loc, :]
    curr_pop_size::Int64 = length(_pop)
    next_pop = @view(ecostate._pop_cache[loc, :])

    # 4 = :sample dimension
    # Using integer indexer for speed (1/2 the allocation and 2/3rds the memory use)
    target_pop_size::Int64 = size(ecostate.pop_sample, 4)

    fill_idx::Int64 = 1

    # Perform a stratified sampling approach, where we keep all members of the largest size
    # class, and then (currently) uniformly sample from the other size classes.
    # Really, we should be estimating a KDE for the lognormal for each size class and then
    # sampling from that

    # Loop from largest to smallest size class
    in_class_idx = falses(length(_pop))
    @floop @inbounds for (sz_id, (sz_end, sz_start)) in enumerate(eachrow(class_diams))
        if sz_id == 1
            # Keep all largest sized corals
            in_class_idx .= sz_start .< _pop
            in_class = count(in_class_idx)
            if in_class == 0
                # No large corals, so skip
                continue
            end

            next_pop[fill_idx:in_class] .= _pop[in_class_idx]
            fill_idx = fill_idx+in_class
            target_pop_size = target_pop_size - in_class
            curr_pop_size = curr_pop_size - in_class
            continue
        end

        # Determine proportion of samples in current size class
        in_class_idx .= sz_start .< _pop .<= sz_end
        in_class = count(in_class_idx)
        if in_class == 0
            continue
        end

        prop_population = in_class / curr_pop_size

        # Note, we take the floor, so we may have less samples than strictly necessary
        # to fill the population sample. We fill these in later.
        n_sample = floor(Int64, target_pop_size * prop_population)
        next_pop[fill_idx:(fill_idx+n_sample-1)] .= sample(
            rng,
            _pop[in_class_idx],
            n_sample
        )
        fill_idx = fill_idx + n_sample
    end

    # Take a weighted sample to fill any remaining empty spots indicated by zero values
    # caused by taking the floor of the population proportion.
    # This should be a small number of the population sample to fill.
    zero_idx::BitVector = next_pop .== 0.0
    if any(zero_idx)
        next_pop[zero_idx] .= rand(rng, _pop, count(zero_idx))

        # Better to take the log-normal within each size class maybe.
        # log_sizes = log.(new_sizes)
        # weights = AverageShiftedHistograms.pdf.([ash(log_sizes; kernel=gaussian)], log_sizes)
        # Below is incorrect - we need to take the probabilities based on interval
        # sampled_indices = sample(
        #     rng,
        #     1:length(new_sizes),
        #     Weights(weights),
        #     count(zero_idx),
        #     replace=false
        # )
    end

    update_sample!(ecostate, ts, loc, grp, next_pop)

    return nothing
end

function generate_example_environment(
    n_timesteps::Int64,
    n_locations::Int64;
    rng::AbstractRNG = Random.GLOBAL_RNG,
    start_year::Int64 = 2020,
    warming_rate::Float32 = 0.05f0,  # °C per year under RCP4.5
    seasonal_amplitude::Float32 = 1.0f0,
    dhw_threshold::Float32 = 2.0f0,  # °C above maximum monthly mean
    noise_amplitude::Float32 = 0.75f0
)
    # Initialize arrays
    dhw_data = zeros(Float32, n_timesteps, n_locations)
    cyclone_data = zeros(Int32, n_timesteps, n_locations)

    # Calculate baseline parameters
    years = range(start_year, length=n_timesteps) .- start_year

    for loc in 1:n_locations
        # Location-specific random offset to create spatial variation
        spatial_offset = rand(rng) * 0.5f0

        # Generate temperature anomaly time series
        for t in 1:n_timesteps
            # Long-term warming trend (RCP4.5-like)
            warming_trend = warming_rate * years[t]

            # Seasonal cycle (simplified)
            seasonal_cycle = seasonal_amplitude * sin(2π * (t % 12) / 12)

            # Random weather variations
            weather_noise = randn(rng) * noise_amplitude

            # Increase variability with warming (supported by climate science)
            variability_scaling = 1.0f0 + 0.1f0 * warming_trend

            # Combined temperature anomaly
            temp_anomaly = (warming_trend +
                          seasonal_cycle * variability_scaling +
                          weather_noise +
                          spatial_offset)

            # Convert temperature anomaly to DHW
            # Only accumulate DHW when anomaly exceeds threshold
            if temp_anomaly > dhw_threshold
                # Accumulate DHW with some decay
                prev_dhw = t > 1 ? dhw_data[t-1, loc] : 0.0f0
                dhw_accumulation = (temp_anomaly - dhw_threshold) / 4.0f0  # Weekly accumulation
                dhw_data[t, loc] = min(12.0f0, prev_dhw * 0.9f0 + dhw_accumulation)
            else
                # DHW decay when below threshold
                prev_dhw = t > 1 ? dhw_data[t-1, loc] : 0.0f0
                dhw_data[t, loc] = max(0.0f0, prev_dhw * 0.7f0)
            end
        end

        # Add occasional extreme events (more frequent in later years)
        n_extreme_events = floor(Int, n_timesteps / 20 * (1.0 + years[end]/50))  # Increasing frequency
        for _ in 1:n_extreme_events
            # Events more likely in later years
            event_time = rand(rng, 1:n_timesteps)
            event_probability = years[event_time] / years[end]

            if rand(rng) < event_probability
                event_duration = rand(rng, 4:8)
                event_magnitude = rand(rng) * 4.0f0 + 4.0f0  # 4-8 DHW

                # Add the extreme event
                for t in event_time:min(event_time + event_duration, n_timesteps)
                    dhw_data[t, loc] = max(dhw_data[t, loc], event_magnitude)
                end
            end
        end
    end
    # dhw_data = zeros(Float32, n_timesteps, n_locations)

    # Generate cyclone data (keeping existing logic)
    # Could be modified to also follow climate projections if needed

    # Combine DHW and cyclone data
    env_data = cat(dhw_data, cyclone_data, dims=3)

    v_axes = (
        Dim{:timestep}(1:n_timesteps),
        Dim{:location}(1:n_locations),
        Dim{:var}([:dhw, :cyclone_category])
    )

    return YAXArray(v_axes, env_data)
end

growth_models = Polynomial[]
model_coefs = Vector{Float32}[]
rmse_scores = Float32[]
r2_scores = Float32[]
for (xi, yi) in zip(eachrow(CoralFlow.bin_edges()), eachrow(CoralFlow.linear_extensions()))
    m = curve_fit(Polynomial, xi, yi, 3)
    push!(rmse_scores, CoralFlow.RMSE(m.(xi), yi))
    push!(r2_scores, CoralFlow.R2(m.(xi), vec(yi)))
    push!(growth_models, m)
    try
        push!(model_coefs, m.coefs)
    catch err
        push!(model_coefs, m.coeffs)
    end
end
