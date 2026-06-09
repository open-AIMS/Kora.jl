import YAXArrays.DD: dims

const SEL = Union{Int64,Colon}  # Defined selector type

"""
    ReefState{F, P, Y3a, Y4a, Y4b}

Mutable container holding the full ecological state of a simulated reef system
across time, space, and functional coral groups.

# Fields
- `wild_population::Array{Vector{F},3}` : Diameter samples (cm) for wild corals,
  indexed `[timestep, location, group]`. Each element is a `Vector` of individual
  colony diameters at that point in the simulation.
- `deployed_population::Array{Vector{F},3}` : Diameter samples (cm) for outplanted
  corals, with the same `[timestep, location, group]` indexing.
- `deployment_times::Array{F,3}` : Number of corals deployed at each
  `[timestep, location, group]` combination.
- `growth_models::Vector{P}` : Per-group growth functions. Each callable maps a
  colony diameter (cm) to its expected diameter at the next annual timestep.
- `survival_models::Vector{P}` : Per-group survival functions. Each callable maps
  a colony diameter (cm) to an annual survival probability.
- `carrying_capacity::Vector{F}` : Maximum coral-bearing area in m^2 for each
  location. Limits total cover and constrains recruitment.
- `depths::Vector{F}` : Water depth in meters at each location. Used to compute
  depth-dependent bleaching mortality coefficients.
- `density::Vector{Int64}` : Maximum colony density in colonies per m^2 for each
  location. Recruitment is suppressed when the total population approaches this
  ceiling.

All other fields (YAXArray tolerance and mortality stores, and fields prefixed
with `_`) are internal implementation details that may change between minor
versions.

# See Also
[`initialize_reef`](@ref), [`run_model!`](@ref)
"""
struct ReefState{
    F<:AbstractFloat,
    P<:Function,
    Y3a<:YAXArray{F,3},
    Y4a<:YAXArray{F,4},
    Y4b<:YAXArray{F,4}
}
    wild_population::Array{Vector{F},3}  # [time, location] ⋅ group
    deployed_population::Array{Vector{F},3}  # [time, location] ⋅ group
    deployment_times::Array{F,3}  # [time, location] ⋅ group
    growth_models::Vector{P}
    survival_models::Vector{P}
    location_scalers::Y3a
    density::Vector{Int64}  # Max population density per location
    depths::Vector{F}       # Depths of each location
    wild_dhw_tolerances::Y4a
    deployed_dhw_tolerances::Y4a
    mortalities::Y4b
    carrying_capacity::Vector{F}  # Area that supports coral in m^2 for each location
    _max_pop_size::Int            # When to sample down population
    _pop_cache::Matrix{F}         # Temporary store of population
    _growth_cache::Matrix{F}      # Temporary store of population growth
    _location_cache::Vector{F}    # Temporary store of population for a given location
    _pop_buffer::Vector{F}        # Oversized buffer for updated population
    _location_buffer::Vector{F}   # Continually updated store of cover for all locations
    _recruit_buffer::Vector{F}    # Continually refreshed store of coral recruits
end

function Base.show(io::IO, ::MIME"text/plain", x::ReefState)
    return println(io, "TODO: Nice printing of ReefState")
end

"""
    Base.copy(rs::ReefState)

Create a thread-safe independent copy of `rs` for use in parallel model
evaluation (e.g. sensitivity analysis with `Threads.@threads`).

`deepcopy(rs)` is unsafe because `YAXArray` is an immutable struct — Julia
returns the same instance rather than creating a new one with copied data.
Fields like `wild_dhw_tolerances` are therefore shared across copies and
corrupted by concurrent writes in `run_model!`.

This method avoids that by:
- Reconstructing each `YAXArray` field explicitly with `copy(field.data)`,
  guaranteeing independent underlying arrays.
- Using `deepcopy` only for `wild_population` / `deployed_population`, whose
  elements are `Vector{F}` objects that must also be independent.
- Using `copy` (shallow) for plain numeric arrays — sufficient since no two
  threads operate on the same reef state instance.
- Sharing `growth_models` and `survival_models` — these hold immutable
  function objects and are safe to read concurrently.
"""
function Base.copy(rs::ReefState)
    return ReefState(
        deepcopy(rs.wild_population),
        deepcopy(rs.deployed_population),
        copy(rs.deployment_times),
        rs.growth_models,
        rs.survival_models,
        YAXArray(dims(rs.location_scalers), copy(rs.location_scalers.data)),
        copy(rs.density),
        copy(rs.depths),
        YAXArray(dims(rs.wild_dhw_tolerances), copy(rs.wild_dhw_tolerances.data)),
        YAXArray(dims(rs.deployed_dhw_tolerances), copy(rs.deployed_dhw_tolerances.data)),
        YAXArray(dims(rs.mortalities), copy(rs.mortalities.data)),
        copy(rs.carrying_capacity),
        rs._max_pop_size,
        copy(rs._pop_cache),
        copy(rs._growth_cache),
        copy(rs._location_cache),
        copy(rs._pop_buffer),
        copy(rs._location_buffer),
        copy(rs._recruit_buffer)
    )
end

@inline function n_timesteps(reef_state::ReefState)::Int64
    return size(reef_state.wild_population, 1)
end

@inline function n_locations(reef_state::ReefState)::Int64
    return size(reef_state.wild_population, 2)
end

@inline function n_groups(reef_state::ReefState)::Int64
    return size(reef_state.wild_population, 3)
end

"""
    total_wild(reef_state::ReefState, ts::Int64, loc::Int64)::Int64

Retrieve the total number of wold corals.
"""
@inline function total_wild(reef_state::ReefState, ts::SEL, loc::SEL)::Int64
    total = 0
    for grp in 1:n_groups(reef_state)
        total += length(reef_state.wild_population[ts, loc, grp])
    end

    return total
end
@inline function total_wild(reef_state::ReefState, ts::SEL, loc::SEL, grp::SEL)::Int64
    return length(reef_state.wild_population[ts, loc, grp])
end

"""
    total_deployed(reef_state::ReefState, ts::Int64, loc::Int64)::Int64

Retrieve the total number of deployed corals.
"""
@inline function total_deployed(reef_state::ReefState, ts::Int64, loc::Int64)::Int64
    total = 0
    for grp in 1:n_groups(reef_state)
        total += total_deployed(reef_state, ts, loc, grp)
    end

    return total
end
@inline function total_deployed(reef_state::ReefState, ts::SEL, loc::SEL, grp::SEL)::Int64
    t::Int64 = try
        length(reef_state.deployed_population[ts, loc, grp])
    catch err
        if !(err isa UndefRefError)
            rethrow(err)
        end

        0
    end

    return t
end

function fill_population_buffer!(
    reef_state::ReefState,
    ts::Int64,
    loc::Int64,
    grp::Int64,
    recruits::Vector{Float32},
    pop_buffer::AbstractVector{Float32}
)::Nothing
    n_wild = total_wild(reef_state, ts, loc, grp)
    pop_buffer[1:n_wild] .= wild_population(reef_state, ts, loc, grp)

    n_deployed = total_deployed(reef_state, ts, loc, grp)
    if n_deployed > 0
        pop_buffer[(n_wild + 1):(n_wild + n_deployed)] .= deployed_population(
            reef_state, ts, loc, grp
        )
    end

    n_recruits = length(recruits)
    if n_recruits > 0
        total_n = n_wild + n_deployed + n_recruits
        if total_n > length(pop_buffer)
            rec_n = length(pop_buffer) - (n_wild + n_deployed)
            pop_buffer[(n_wild + n_deployed + 1):end] .= recruits[1:rec_n]
        else
            pop_buffer[(n_wild + n_deployed + 1):(n_wild + n_deployed + n_recruits)] .=
                recruits
        end
    end

    return nothing
end

"""
    total_population(reef_state::ReefState, ts::Int64, loc::Int64)::Int64

Retrieve the total population count.
"""
@inline function total_population(reef_state::ReefState, ts::Int64, loc::Int64)::Int64
    total = 0
    for grp in 1:n_groups(reef_state)
        total += length(reef_state.wild_population[ts, loc, grp])
        total += length(reef_state.deployed_population[ts, loc, grp])
    end

    return total
end
@inline function total_population(reef_state::ReefState, ts::SEL, loc::SEL, grp::SEL)::Int64
    wild = length(reef_state.wild_population[ts, loc, grp])
    deployed = length(reef_state.deployed_population[ts, loc, grp])
    return wild + deployed
end

@inline function wild_population(
    reef_state::ReefState, ts::SEL, loc::SEL, grp::SEL
)::SubArray
    return @inbounds @views reef_state.wild_population[ts, loc, grp][:]
end

@inline function deployed_population(
    reef_state::ReefState, ts::SEL, loc::SEL, grp::SEL
)::SubArray
    return @inbounds @views reef_state.deployed_population[ts, loc, grp][:]
end

@inline function coral_population(
    reef_state::ReefState, ts::SEL, loc::SEL, grp::SEL;
    cache=@view(reef_state._pop_cache[loc, :])
)::SubArray
    return coral_population!(reef_state, ts, loc, grp, cache)
end
function coral_population!(
    reef_state::ReefState, ts::SEL, loc::SEL, grp::SEL, cache::AbstractVector
)::SubArray
    # TODO: Handle request for sets of groups or locations
    n_wild = total_wild(reef_state, ts, loc, grp)
    n_deployed = total_deployed(reef_state, ts, loc, grp)

    cache[1:n_wild] .= wild_population(reef_state, ts, loc, grp)

    if n_deployed > 0
        cache[(n_wild + 1):(n_wild + n_deployed)] .= deployed_population(
            reef_state, ts, loc, grp
        )
    end

    return @inbounds @view(cache[1:(n_wild + n_deployed)])
end

function update_wild_sample!(
    reef_state::ReefState, ts::SEL, loc::SEL, group::SEL, pop::AbstractVector
)::Nothing
    @inbounds reef_state.wild_population[ts, loc, group] = pop

    return nothing
end

function update_deployed_sample!(
    reef_state::ReefState, ts::SEL, loc::SEL, group::SEL, pop::AbstractVector
)::Nothing
    @inbounds reef_state.deployed_population[ts, loc, group] = pop

    return nothing
end

function update_dhw_tol_mean!(
    reef_state::ReefState, ts::Int64, grp::Int64, vals::AbstractVector{Float32}
)::Nothing
    reef_state.wild_dhw_tolerances.data[ts, :, grp, 1] = vals

    return nothing
end
function update_dhw_tol_mean!(
    reef_state::ReefState, ts::Int64, loc::Int64, grp::Int64, val::Float32
)::Nothing
    reef_state.wild_dhw_tolerances.data[ts, loc, grp, 1] = val

    return nothing
end

function update_dhw_tol_std!(
    reef_state::ReefState, ts::Int64, loc::Int64, grp::Int64, val::Float32
)::Nothing
    reef_state.wild_dhw_tolerances.data[ts, loc, grp, 2] = val

    return nothing
end

"""
    initialize_reef(;
        n_timesteps::Int=75,
        n_locs::Int=100,
        group_names=Kora.TARGET_GROUPS,
        density::Union{Int64,Vector{Int64}}=20,
        area=90.0,
        depths::Union{Float64,Vector{Float64}}=9.0,
        growth_models::AbstractCoralBehavior=Kora.growth_models,
        survival_models::AbstractCoralBehavior=Kora.survival_models
    )::ReefState

Allocate and return a `ReefState` sized for `n_locs` reef locations and
`n_timesteps` annual time steps.

All population arrays are initialised empty. Call `initialize_coral_population!`
to seed the starting population before running a simulation.

# Arguments
- `n_timesteps` : Number of annual time steps to allocate (default: `75`).
- `n_locs` : Number of reef locations (default: `100`).
- `group_names` : Labels for the functional coral groups. Must match the groups
  used to fit `growth_models` and `survival_models`
  (default: `Kora.TARGET_GROUPS`, the five groups used by the bundled models).
- `density` : Maximum colony density in colonies per m^2. Provide a scalar to
  apply the same ceiling to every location, or a per-location `Vector{Int64}`
  (default: `20`).
- `area` : Reef area available for coral cover in m^2. Provide a scalar or a
  per-location `Vector` (default: `90.0`).
- `depths` : Water depth in meters. Provide a scalar or a per-location
  `Vector{Float64}`. Depth controls bleaching mortality coefficients
  (default: `9.0`).
- `growth_models` : Fitted growth model collection, one function per functional
  group. Defaults to the package-level offshore-north models loaded from the
  bundled JSON asset at package load time.
- `survival_models` : Fitted survival model collection, one function per
  functional group. Defaults to the package-level offshore-north models loaded
  from the bundled JSON asset at package load time.

# Returns
`ReefState` : An empty reef state ready for population initialisation.

# Examples
```jldoctest
julia> using Kora

julia> rs = initialize_reef(; n_timesteps=10, n_locs=3);

julia> n_timesteps(rs), n_locations(rs)
(10, 3)
```

# See Also
[`initialize_coral_population!`](@ref), [`run_model!`](@ref),
[`load_models`](@ref)
"""
function initialize_reef(;
    n_timesteps=75,
    n_locs=100,
    group_names=Kora.TARGET_GROUPS,
    density::Union{Int64,Vector{Int64}}=20,  # Max density per unit area
    area=90.0,
    depths::Union{Float64,Vector{Float64}}=9.0,
    growth_models::AbstractCoralBehavior=Kora.growth_models,
    survival_models::AbstractCoralBehavior=Kora.survival_models
)
    n_groups = length(group_names)

    loc_ax = (
        Dim{:scaler}([:growth, :mortality]),
        Dim{:location}(1:n_locs),
        Dim{:group}(group_names)
    )
    location_scalers = YAXArray(loc_ax, fill(1.0f0, 2, n_locs, n_groups))

    # Initialize empty vector arrays
    wild_population = Array{Vector{Float32},3}(undef, n_timesteps, n_locs, n_groups)
    deployed_population = Array{Vector{Float32},3}(undef, n_timesteps, n_locs, n_groups)

    wild_population .= [Float32[]]
    deployed_population .= [Float32[]]
    deployment_times = zeros(Float32, n_timesteps, n_locs, n_groups)

    dhw_ax = (
        Dim{:timestep}(1:n_timesteps),
        Dim{:location}(1:n_locs),
        Dim{:group}(group_names),
        Dim{:factor}([:mean, :stdev])
    )
    wild_dhw_tol = YAXArray(dhw_ax, zeros(Float32, n_timesteps, n_locs, n_groups, 2))
    deployed_dhw_tol = YAXArray(dhw_ax, zeros(Float32, n_timesteps, n_locs, n_groups, 2))

    mort_names = [:dhw, :cyclone]
    mort_ax = (
        Dim{:timestep}(1:n_timesteps),
        Dim{:location}(1:n_locs),
        Dim{:group}(group_names),
        Dim{:mortality}(mort_names)
    )
    mort = YAXArray(
        mort_ax, zeros(Float32, n_timesteps, n_locs, n_groups, length(mort_names))
    )

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
        throw(
            ArgumentError(
                "Area must be a single value or a vector indicating area of each location"
            )
        )
    end

    if density isa Integer
        density = fill(density, n_locs)
    elseif !(area isa Vector{<:Real})
        msg = "Density must be a single value or a vector indicating maximum population density of each location"
        throw(ArgumentError(msg))
    end

    buffer_size = Int64(maximum(area) * maximum(density) * 5)
    reef_state = ReefState(
        wild_population,
        deployed_population,
        deployment_times,
        growth_models.models,
        survival_models.models,
        location_scalers,
        density,
        Float32.(depths),
        wild_dhw_tol,
        deployed_dhw_tol,
        mort,
        area,
        10_000,
        zeros(Float32, n_locs, buffer_size),
        zeros(Float32, n_locs, buffer_size),
        zeros(Float32, n_locs),
        zeros(Float32, buffer_size),
        zeros(Float32, n_locs),
        zeros(Float32, n_locs)
    )

    return reef_state
end

"""
    reset!(reef_state::ReefState)

Clear out existing results with empty vector arrays, keeping only the initial state.
"""
function reset!(reef_state::ReefState)
    reef_state.wild_population[2:end, :, :] .= [Float32[]]
    reef_state.deployed_population[2:end, :, :] .= [Float32[]]

    reef_state.wild_dhw_tolerances[2:end, :, :, 1] .= 0.0f0
    return reef_state.deployed_dhw_tolerances[2:end, :, :, 1] .= 0.0f0
end

function size_distribution()::Vector{LogNormal}
    return [
        LogNormal{Float32}(2.2382145f0, 0.74870664f0),  # tabular Acropora
        LogNormal{Float32}(2.1610258f0, 0.64033926f0),  # corymbose Acropora
        LogNormal{Float32}(1.7919482f0, 0.62202924f0),  # Pocillopora + non-Acropora corymbose
        LogNormal{Float32}(2.0454452f0, 1.4701339f0),   # Small massives and encrusting
        LogNormal{Float32}(1.9826132f0, 1.4488107f0)    # Large massives
    ]
end

"""
    initialize_coral_population!(
        reef_state::ReefState,
        loc::Int64,
        target_pop_size::Int64;
        group_proportions::Vector{Float32}=[0.1f0, 0.2f0, 0.25f0, 0.2f0, 0.25f0],
        rng::AbstractRNG=Random.GLOBAL_RNG
    )::Nothing

    initialize_coral_population!(
        reef_state::ReefState;
        rng::AbstractRNG=Random.GLOBAL_RNG
    )

Seed the initial coral population at timestep 1 with colony diameter samples
drawn from per-group log-normal size distributions.

The no-location convenience method seeds all locations using a target population
size derived from the maximum carrying capacity in `reef_state`. Colony sizes are
drawn from group-specific truncated log-normal distributions and written as
diameter vectors (cm) into `wild_population[1, loc, grp]`.

# Arguments
- `reef_state` : The `ReefState` to populate. Population data are written
  in-place at timestep 1.
- `loc` : Index of the location to seed (1-based).
- `target_pop_size` : Total number of colonies to place at this location. Actual
  per-group counts are `round(target_pop_size * group_proportions[grp])`.
- `group_proportions` : Proportion of the total population assigned to each
  functional group. Must sum to 1.0 (checked with `atol=1e-6`).
  Default: `[0.10, 0.20, 0.25, 0.20, 0.25]` matching the five bundled groups.
- `rng` : Random number generator for reproducible diameter draws
  (default: `Random.GLOBAL_RNG`).

# Returns
`Nothing`

# See Also
[`initialize_reef`](@ref), [`run_model!`](@ref)
"""
function initialize_coral_population!(
    reef_state::ReefState,
    loc::Int64,
    target_pop_size::Int64;
    group_proportions::Vector{Float32}=[0.1f0, 0.2f0, 0.25f0, 0.2f0, 0.25f0],
    size_dist=size_distribution(),
    rng::AbstractRNG=Random.GLOBAL_RNG
)::Nothing
    # Verify proportions sum to 1
    if !isapprox(sum(group_proportions), 1.0f0; atol=1e-6)
        throw(ArgumentError("Group proportions must sum to 1.0"))
    end

    edges = bin_edges()
    for grp in 1:n_groups(reef_state)
        dist = truncated(size_dist[grp], 0.0, maximum(edges[grp, :]))

        # Scale the number of samples by the desired proportion
        n_samples = round(Int, target_pop_size * group_proportions[grp])
        initial_population = convert.(Float32, rand(rng, dist, n_samples))

        update_wild_sample!(reef_state, 1, loc, grp, initial_population)
    end

    reef_state.wild_dhw_tolerances[1, :, 1, At(:mean)] .= 3.751612251  # tabular Acropora
    reef_state.wild_dhw_tolerances[1, :, 2, At(:mean)] .= 4.081622683  # corymbose Acropora
    reef_state.wild_dhw_tolerances[1, :, 3, At(:mean)] .= 4.487465256  # Pocillopora + non-Acropora corymbose
    reef_state.wild_dhw_tolerances[1, :, 4, At(:mean)] .= 6.165751937  # Small massives and encrusting
    reef_state.wild_dhw_tolerances[1, :, 5, At(:mean)] .= 7.153507902  # Large massives

    reef_state.wild_dhw_tolerances[:, :, 1, At(:stdev)] .= 2.904433676  # tabular Acropora
    reef_state.wild_dhw_tolerances[:, :, 2, At(:stdev)] .= 3.159922076  # corymbose Acropora
    reef_state.wild_dhw_tolerances[:, :, 3, At(:stdev)] .= 3.474118416  # Pocillopora + non-Acropora corymbose
    reef_state.wild_dhw_tolerances[:, :, 4, At(:stdev)] .= 4.773419097  # Small massives and encrusting
    reef_state.wild_dhw_tolerances[:, :, 5, At(:stdev)] .= 5.538122776  # Large massives

    return nothing
end
function initialize_coral_population!(
    reef_state::ReefState;
    rng::AbstractRNG=Random.GLOBAL_RNG
)
    n_locs = n_locations(reef_state)
    sample_size = ceil(Int64, maximum(reef_state.carrying_capacity) * 5)
    for loc in 1:n_locs
        initialize_coral_population!(reef_state, loc, sample_size; rng=rng)
    end
end

"""
    deploy_corals!(reef_state, ts, loc, n, grp; rng=Random.GLOBAL_RNG)

Seed `n` outplanted coral colonies of functional group `grp` at location `loc`
and timestep `ts`.

Colony initial diameters are sampled from the truncated log-normal size distribution
for `grp` (bounded by the group's diameter bin edges) and stored in
`reef_state.deployed_population[ts, loc, grp]`. Any existing deployed population
at that slot is overwritten.

This function is called internally when coral deployment is active in `run_model!`.
Direct calls are supported for testing but are not part of the standard simulation
workflow -- use the deployment configuration in `run_model!` to trigger outplanting
within a simulation run.

# Arguments
- `reef_state` : `ReefState` to update in-place.
- `ts` : Timestep index at which deployment occurs (1-based).
- `loc` : Location index at which corals are deployed (1-based).
- `n` : Number of colonies to deploy.
- `grp` : Functional group index (1-based).
- `rng` : Random number generator for reproducible diameter draws
  (default: `Random.GLOBAL_RNG`).

# Returns
`Nothing`

# See Also
[`initialize_coral_population!`](@ref), [`run_model!`](@ref)
"""
function deploy_corals!(reef_state, ts, loc, n, grp; rng=Random.GLOBAL_RNG)
    size_dist = size_distribution()[grp]
    edges = bin_edges()[grp, :]

    dist = truncated(size_dist, 0.0, maximum(edges[grp, :]))
    deploy_sample = convert.(Float32, rand(rng, dist, n))

    return update_deployed_sample!(reef_state, ts, loc, grp, deploy_sample)
end

function update_pop_cache!(reef_state, current_diams, loc)::Nothing
    next_pop = @view(reef_state._pop_cache[loc, :])

    reef_state._pop_cache[loc, 1:length(current_diams)] .= current_diams
    if length(current_diams) < length(next_pop)
        next_pop[(length(current_diams) + 1):end] .= 0.0f0
    end

    return nothing
end

"""
    resample_wild_population!(
        reef_state::ReefState,
        ts::Int64,
        loc::Int64,
        grp::Int64,
        diams::AbstractArray{Float32},
        class_diams::Matrix{Float32};
        rng::AbstractRNG=Random.GLOBAL_RNG
    )::Nothing

Sample from population to avoid exceeding density threshold (coral per m²).
"""
function resample_wild_population!(
    reef_state::ReefState,
    ts::Int64,
    loc::Int64,
    grp::Int64,
    diams::AbstractArray{Float32},
    class_diams::Matrix{Float32};
    rng::AbstractRNG=Random.GLOBAL_RNG
)::Nothing
    next_pop = @view(reef_state._pop_cache[loc, :])
    max_sample_size =
        Int64(reef_state.carrying_capacity[loc] * reef_state.density[loc]) * 20

    # If population is below max sample size, keep everything
    if length(diams) <= max_sample_size
        update_pop_cache!(reef_state, diams)
    else
        # Only resample if we exceed max sample size
        # Use size-stratified sampling to maintain size distribution
        size_weights = zeros(Float32, length(diams))
        for (sz_end, sz_start) in eachrow(class_diams)
            in_class = sz_start .< diams .<= sz_end
            if any(in_class)
                # Weight by colony size to maintain proper cover representation
                size_weights[in_class] .= diams[in_class] ./ sum(diams[in_class])
            end
        end

        # Sample proportionally to maintain size structure
        sampled_idx = sample(
            rng, 1:length(diams), Weights(size_weights), max_sample_size; replace=false
        )
        next_pop .= diams[sampled_idx]
    end

    update_wild_sample!(reef_state, ts, loc, grp, next_pop[next_pop .> 0.0])

    return nothing
end

"""
    coral_cover(reef_state::ReefState)::Vector{Float32}
    coral_cover(reef_state::ReefState, ts::Int64)::Vector{Float32}
    coral_cover(reef_state::ReefState, ts::Int64, loc::Int64)::Float32
    coral_cover(diams::AbstractVector{<:AbstractFloat})

Compute total coral cover in m^2, summing over all wild and deployed colonies
across all functional groups.

The no-timestep form returns a `Vector` of summed cover across all locations for
each timestep. The single-timestep form returns a per-location `Vector` at `ts`.
The two-argument form returns a scalar for one timestep and one location. The
bare-vector form computes cover from a diameter vector (cm) directly, without
requiring a `ReefState`.

Colony area is computed as pi/4 * (d/100)^2 (m^2) for each diameter d in cm.

# Arguments
- `reef_state` : Source of population data.
- `ts` : Timestep index (1-based).
- `loc` : Location index (1-based).
- `diams` : Vector of colony diameters in cm.

# Returns
`Vector{Float32}` or `Float32` : Coral cover in m^2.

# See Also
[`group_cover`](@ref), [`juvenile_cover`](@ref), [`cover_cm_to_m2`](@ref)
"""
function coral_cover(reef_state::ReefState)::Vector{Float32}
    covers = zeros(Float32, n_timesteps(reef_state))
    for ts in axes(covers, 1)
        @inbounds covers[ts] = sum(coral_cover(reef_state, ts))
    end

    return covers
end
function coral_cover(reef_state::ReefState, ts::Int64)::Vector{Float32}
    loc_covers = reef_state._location_buffer
    for loc in eachindex(loc_covers)
        @inbounds loc_covers[loc] = coral_cover(reef_state, ts, loc)
    end

    return loc_covers
end
function coral_cover(reef_state::ReefState, ts::Int64, loc::Int64)::Float32
    total = 0.0f0
    for grp in 1:n_groups(reef_state)
        total += sum(cover_cm_to_m2(reef_state.wild_population[ts, loc, grp]))
        total += sum(cover_cm_to_m2(reef_state.deployed_population[ts, loc, grp]))
    end

    return total
end
function coral_cover(diams::AbstractVector{<:AbstractFloat})
    return sum(cover_cm_to_m2.(diams))
end

"""
    group_cover(reef_state::ReefState)::Matrix{Float32}
    group_cover(reef_state::ReefState, ts::Int64)::Vector{Float32}

Compute mean coral cover in m^2 per functional group, averaged across all
locations.

The no-timestep form delegates to `group_cover_timeseries` and returns results
for every timestep as a matrix. The single-timestep form returns a per-group
vector at `ts`.

# Arguments
- `reef_state` : Source of population data.
- `ts` : Timestep index (1-based).

# Returns
`Matrix{Float32}` with shape `(n_timesteps, n_groups)`, or `Vector{Float32}`
with one element per functional group. Values are mean cover in m^2 across all
locations.

# See Also
[`coral_cover`](@ref), [`group_cover_timeseries`](@ref),
[`juvenile_cover`](@ref)
"""
function group_cover(reef_state::ReefState)::Matrix{Float32}
    return group_cover_timeseries(reef_state)
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

"""
    group_cover_timeseries(reef_state::ReefState)::Matrix{Float32}

Compute mean coral cover in m^2 per functional group for every timestep,
averaged across all locations.

# Arguments
- `reef_state` : Source of population data.

# Returns
`Matrix{Float32}` : Cover matrix of shape `(n_timesteps, n_groups)`. Rows are
timesteps; columns are functional groups in the order defined by
`reef_state.wild_population`.

# See Also
[`group_cover`](@ref), [`coral_cover`](@ref)
"""
function group_cover_timeseries(reef_state::ReefState)::Matrix{Float32}
    n_ts = n_timesteps(reef_state)
    n_grp = n_groups(reef_state)
    covers = zeros(Float32, n_ts, n_grp)

    for ts in 1:n_ts
        covers[ts, :] = group_cover(reef_state, ts)
    end

    return covers
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
    for grp in 1:n_grps, loc in 1:n_locs
        # Get population sample for this location/group
        pop = coral_population(reef_state, ts, loc, grp)

        # Sum cover of colonies above maturity threshold
        mature_cover[loc, grp] = cover_cm_to_m2(pop[pop .>= thresholds[grp]])
    end

    return mature_cover
end

"""
    recruit_cover(recruits::Array{Float32})::Vector{Float32}
    recruit_cover(recruits::Matrix{Vector{Float32}})::Vector{Float32}
    recruit_cover(ecostate::ReefState, recruits::Array{Float32})
    recruit_cover(ecostate::ReefState, recruits::Matrix{Vector{Float32}})

Compute total coral cover in m^2 for a cohort of new recruits, summed per location.

The single-argument forms allocate a fresh output vector. The two-argument forms
write into `ecostate._recruit_buffer` and return that view, avoiding allocation
inside the simulation loop.

Colony area is computed as `pi/4 * (d/100)^2` (m^2) for each recruit diameter `d`
in cm via [`cover_cm_to_m2`](@ref).

# Arguments
- `ecostate` : `ReefState` whose `_recruit_buffer` is used to store results.
- `recruits` : Per-location recruit diameter data. Either a 3-D `Array{Float32}`
  with axes `[location, group, colony]`, or a `Matrix{Vector{Float32}}` with
  dimensions `[location, group]`.

# Returns
`Vector{Float32}` : Total recruit cover in m^2, one element per location.

# See Also
[`coral_cover`](@ref), [`cover_cm_to_m2`](@ref)
"""
function recruit_cover(recruits::Array{Float32})::Vector{Float32}
    loc_cover = zeros(Float32, axes(recruits, 1))  # zeros for each location
    tmp_cover = cover_cm_to_m2.(recruits)
    @inbounds for i in eachindex(loc_cover)
        loc_cover[i] = sum(@view(tmp_cover[i, :, :]))
    end

    return loc_cover
end
function recruit_cover(recruits::Matrix{Vector{Float32}})::Vector{Float32}
    loc_cover = zeros(Float32, axes(recruits, 1))  # zeros for each location
    tmp_cover = cover_cm_to_m2.(recruits)
    @inbounds for i in eachindex(loc_cover)
        loc_cover[i] = sum(@view(tmp_cover[i, :, :]))
    end

    return loc_cover
end

function recruit_cover(ecostate::ReefState, recruits::Array{Float32})
    loc_cover = ecostate._recruit_buffer

    tmp_cover = cover_cm_to_m2.(recruits)
    @inbounds for i in eachindex(loc_cover)
        loc_cover[i] = sum(@view(tmp_cover[i, :, :]))
    end

    return loc_cover
end
function recruit_cover(ecostate::ReefState, recruits::Matrix{Vector{Float32}})
    loc_cover = ecostate._recruit_buffer

    tmp_cover = cover_cm_to_m2.(recruits)
    @inbounds for i in eachindex(loc_cover)
        loc_cover[i] = sum(@view(tmp_cover[i, :]))
    end

    return loc_cover
end

# Non-bootstrapped juvenile cover functions
"""
    juvenile_cover(
        reef_state::ReefState,
        ts::Int64;
        juvenile_threshold::Union{Nothing,Float32,Vector{Float32}}=nothing
    )::Vector{Float32}

Compute mean cover of sub-mature corals in m^2 per functional group at a single
timestep, averaged across all locations.

A colony is classified as juvenile when its diameter is strictly less than the
maturity threshold for its group. The default thresholds come from
`Kora.mature_size_thresholds()`; pass a custom value to override.

# Arguments
- `reef_state` : Source of population data.
- `ts` : Timestep index (1-based).
- `juvenile_threshold` : Diameter threshold in cm below which a colony counts
  as juvenile. Provide a scalar `Float32` to apply the same value to all groups,
  a `Vector{Float32}` for per-group thresholds, or `nothing` to use the
  package defaults (default: `nothing`).

# Returns
`Vector{Float32}` : Mean juvenile cover in m^2, one element per functional
group.

# See Also
[`group_cover`](@ref), [`coral_cover`](@ref)
"""
function juvenile_cover(
    reef_state::ReefState,
    ts::Int64;
    juvenile_threshold::Union{Nothing,Float32,Vector{Float32}}=nothing
)::Vector{Float32}
    n_grp = n_groups(reef_state)
    n_loc = n_locations(reef_state)
    means = zeros(Float32, n_grp)

    if isnothing(juvenile_threshold)
        mature_sizes = mature_size_thresholds()
    else
        if !(juvenile_threshold isa Vector)
            mature_sizes = fill(juvenile_threshold, n_grp)
        else
            mature_sizes = juvenile_threshold
        end
    end

    for grp in 1:n_grp
        total_cover = 0.0f0
        grp_mature = mature_sizes[grp]
        for loc in 1:n_loc
            pop = coral_population(reef_state, ts, loc, grp)
            juveniles = pop[0.0f0 .< pop .< grp_mature]
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
    juvenile_threshold::Union{Nothing,Float32,Vector{Float32}}=nothing
)::Matrix{Float32}
    n_ts = n_timesteps(reef_state)
    n_grp = n_groups(reef_state)
    means = zeros(Float32, n_ts, n_grp)

    if isnothing(juvenile_threshold)
        mature_sizes = mature_size_thresholds()
    else
        if !(juvenile_threshold isa Vector)
            mature_sizes = fill(juvenile_threshold, n_grp)
        else
            mature_sizes = juvenile_threshold
        end
    end

    for ts in 1:n_ts
        means[ts, :] = juvenile_cover(reef_state, ts; juvenile_threshold=mature_sizes)
    end

    return means
end

"""
    generate_example_environment(
        n_years::Int64,
        n_locations::Int64;
        rng::AbstractRNG=Random.GLOBAL_RNG,
        start_year::Int64=2020,
        with_dhw=true,
        warming_rate::Float32=0.15f0,
        seasonal_amplitude::Float32=1.2f0,
        dhw_threshold::Float32=4.0f0,
        noise_amplitude::Float32=0.9f0
    )

Generate synthetic environmental data for coral reef modeling, specifically Degree Heating
Weeks (DHW) trajectories that simulate realistic marine heatwave patterns under climate
change scenarios.

# Extended help

Generate plausible DHW time series by combining multiple environmental components:

1. **Long-term warming trend**: Simulates gradual ocean warming (like RCP4.5 scenarios)
2. **Seasonal cycles**: Models natural temperature variations throughout the year
3. **Weather noise**: Adds realistic short-term temperature fluctuations
4. **Spatial variation**: Creates location-specific temperature offsets
5. **Acute heatwave events**: Superimposes extreme marine heatwave events
6. **DHW accumulation**: Converts temperature anomalies to ecologically-relevant DHW values

The DHW calculation follows coral bleaching research where:
- DHW accumulates when temperatures exceed a threshold (default 4°C above baseline)
- Values decay over time when temperatures drop
- Extreme events can cause rapid DHW spikes that lead to mass bleaching

## Temperature Anomaly Construction
For each location and timestep, the temperature anomaly is built from:

```julia
temp_anomaly = warming_trend + seasonal_cycle + weather_noise + spatial_offset
```

Where:
- **warming_trend**: Linear increase over time (`warming_rate * years_elapsed`)
- **seasonal_cycle**: Sinusoidal pattern with amplitude that increases over time
- **weather_noise**: Random normal variations that get larger in later years
- **spatial_offset**: Location-specific random offset (0 - 0.8°C)

## DHW Accumulation Rules
- **Above threshold**: DHW accumulates as `(temp_anomaly - threshold) / 4.0`
- **Below threshold**: DHW decays rapidly (`previous_DHW * 0.7`)
- **During heating**: DHW decays slowly (`previous_DHW * 0.92`)
- **Soft cap**: Values above 20 DHW are dampened but can still fluctuate

## Acute Event Generation
The function adds realistic extreme heatwave events:
- **Frequency**: ~1 event per year on average (`n_timesteps / 12`)
- **Probability**: Increases over time (simulating worsening climate)
- **Duration**: 2-5+ weeks, longer in later years
- **Intensity**: 8-25+ DHW, stronger in later years
- **Shape**: Rapid onset (30% of duration) → peak → gradual decline (30% of duration)

## Ecological Realism
The parameters are tuned based on coral bleaching research:
- **4 DHW**: Threshold where bleaching typically begins
- **8+ DHW**: Significant bleaching and mortality expected
- **20+ DHW**: Severe bleaching events (soft-capped with fluctuations)

# Arguments
- `n_years`: Number of time steps to simulate (in years)
- `n_locations`: Number of spatial locations across the reef system
- `rng`: Random number generator for reproducible results
- `start_year`: Starting year for the simulation (default: 2020)
- `with_dhw`: Whether to generate DHW data or return zeros (default: true)
- `warming_rate`: Rate of long-term warming per year in °C (default: 0.15°C/year)
- `seasonal_amplitude`: Strength of seasonal temperature cycles (default: 1.2°C)
- `dhw_threshold`: Temperature threshold above which DHW accumulates (default: 4.0°C)
- `noise_amplitude`: Magnitude of random weather variations (default: 0.9°C)
"""
function generate_example_environment(
    n_years::Int64,
    n_locations::Int64;
    rng::AbstractRNG=Random.GLOBAL_RNG,
    start_year::Int64=2020,
    with_dhw=true,
    warming_rate::Float32=0.15f0,
    seasonal_amplitude::Float32=1.2f0,
    dhw_threshold::Float32=4.0f0,
    noise_amplitude::Float32=0.9f0
)::YAXArray
    # Calculate baseline parameters
    years = range(start_year; length=n_years) .- start_year

    # Initialize vectors
    dhw_data = zeros(Float32, n_years, n_locations)

    if with_dhw
        for loc in 1:n_locations
            # Location-specific random offset to create spatial variation
            spatial_offset = rand(rng) * 0.8f0  # Increased spatial variation

            # Track consecutive weeks above threshold
            weeks_above_threshold = 0

            # Generate temperature anomaly time series
            for t in 1:n_years
                # Long-term warming trend
                warming_trend = warming_rate * years[t]

                # Seasonal cycle with increasing amplitude over time
                seasonal_scaling = 1.0f0 + 0.2f0 * (years[t] / years[end])
                seasonal_cycle =
                    seasonal_amplitude * seasonal_scaling * sin(2π * (t % 12) / 12)

                # Random weather variations with increased variability over time
                weather_scaling = 1.0f0 + 0.3f0 * (years[t] / years[end])
                weather_noise = randn(rng) * noise_amplitude * weather_scaling

                # Combined temperature anomaly
                temp_anomaly = (
                    warming_trend +
                    seasonal_cycle +
                    weather_noise +
                    spatial_offset
                )

                prev_dhw = t > 1 ? dhw_data[t - 1, loc] : 0.0f0

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
                    # Reduce value for slower decay during heating
                    base_dhw = prev_dhw * 0.92f0
                    raw_dhw = base_dhw + dhw_accumulation

                    # Apply soft cap with fluctuations
                    soft_cap = 20.0f0
                    if raw_dhw > soft_cap
                        excess = raw_dhw - soft_cap
                        damping_factor = 1.0f0 / (1.0f0 + 0.3f0 * excess)

                        # Larger fluctuations in later years
                        time_factor = years[t] / years[end]
                        max_fluctuation = 4.0f0 + 2.0f0 * time_factor
                        fluctuation =
                            (rand(rng) - 0.5f0) * min(max_fluctuation, excess * 0.6f0)

                        dhw_data[t, loc] =
                            soft_cap + (excess * damping_factor) + fluctuation
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
            # lower denominator for more events
            n_extreme_events = floor(Int, n_years / 12)

            for _ in 1:n_extreme_events
                event_time = rand(rng, 1:n_years)
                time_progress = years[event_time] / years[end]
                event_probability = time_progress * 1.8f0  # Higher probability in later years

                if rand(rng, Float32) < event_probability
                    # Duration increases with time
                    base_duration = 2:5
                    extra_duration = rand(rng, 0:floor(Int, 3 * time_progress))
                    event_duration = rand(rng, base_duration) + extra_duration

                    # Base magnitude increases with time (8 - 17 DHW)
                    base_magnitude = 8.0f0 + 17.0f0 * time_progress
                    event_magnitude = rand(rng) * 5.0f0 + base_magnitude

                    # Add the extreme event with realistic onset/decline
                    for t in event_time:min(event_time + event_duration, n_years)
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
    end

    # Variable axes
    v_axes = (
        Dim{:timestep}(1:n_years),
        Dim{:location}(1:n_locations),
        Dim{:variable}([:dhw])
    )

    return YAXArray(v_axes, reshape(dhw_data, n_years, n_locations, 1))
end

"""
    generate_environment(dhw::Matrix{Float32}; start_year::Int64=2020)::YAXArray
    generate_environment(dhw::YAXArray; start_year::Int64=2020)::YAXArray

Wrap user-supplied Degree Heating Weeks data in a correctly structured environment
YAXArray, suitable for direct use with Kora.jl model runs.

`n_years` and `n_locations` are inferred from the input dimensions. Structure
creation is delegated to `generate_example_environment` so dimension names and
axis labels are always consistent. The synthetic DHW values produced by that
function are then replaced with the caller's real data.

# Arguments
- `dhw` : DHW data with shape `(n_timesteps, n_locations)`. Each row is one
  year (or timestep) and each column is one reef location. Accepts either a
  `Matrix{Float32}` or a 2D `YAXArray` whose first dimension is timestep and
  second dimension is location.
- `start_year` : First year label for the timestep axis (default: 2020). Passed
  through to `generate_example_environment` for axis labelling only; it does not
  alter the data.

# Returns
A 3D `YAXArray` with axes `(Dim{:timestep}, Dim{:location}, Dim{:variable})`
identical in structure to the output of `generate_example_environment`, with the
`:dhw` variable populated from `dhw`.

# Notes
Two advisory warnings are issued (not errors). A minimum-floor check fires when
`minimum(dhw) > 20`, because real DHW data always contains near-zero values
during non-bleaching periods -- a uniformly high floor is the primary indicator
that raw sea-surface temperature (~25-32 degrees C) was passed instead of DHW.
A ceiling check fires when `maximum(dhw) > 40`, approximately twice the ~20 DHW
projected under SSP5-8.5; values above this threshold are likely a data quality
issue rather than an intentional scenario.

# Examples
```jldoctest
julia> using Kora

julia> dhw = zeros(Float32, 10, 5);

julia> env = generate_environment(dhw);

julia> size(env)
(10, 5, 1)
```

# See Also
[`generate_example_environment`](@ref)
"""
function generate_environment(dhw::Matrix{Float32}; start_year::Int64=2020)::YAXArray
    n_years, n_locs = size(dhw)

    if n_years == 0 || n_locs == 0
        throw(ArgumentError(
            "DHW matrix must have at least one timestep and one location, " *
            "got size $(size(dhw))"
        ))
    end

    if any(dhw .< 0.0f0)
        throw(ArgumentError(
            "DHW values must be non-negative; found values below zero"
        ))
    end

    # A high minimum across all locations and timesteps is the primary indicator
    # that raw SST (typically 25-32 degrees C, never near zero) has been supplied
    # instead of DHW (which should contain many near-zero values during non-bleaching
    # periods and early simulation years).
    if minimum(dhw) > 20.0f0
        @warn "Minimum DHW value is $(minimum(dhw)) across all locations and timesteps. " *
              "Real DHW data should include near-zero values during non-bleaching periods. " *
              "A uniformly high floor may indicate raw sea-surface temperature was passed " *
              "instead of accumulated degree heating weeks."
    end

    # Values above 40 DHW are approximately twice the ~20 DHW projected under SSP5-8.5
    # and are likely a data quality issue rather than an intentional scenario.
    if maximum(dhw) > 40.0f0
        @warn "DHW values exceed 40 (maximum = $(maximum(dhw))). This is roughly twice the " *
              "~20 DHW projected under SSP5-8.5. Consider checking for data quality issues " *
              "or confirming the scenario is intentionally extreme."
    end

    env = generate_example_environment(
        n_years, n_locs; start_year=start_year, with_dhw=false
    )

    env.data[:, :, 1] .= dhw

    return env
end

function generate_environment(dhw::YAXArray; start_year::Int64=2020)::YAXArray
    if ndims(dhw) != 2
        throw(ArgumentError(
            "DHW YAXArray must be 2D with dimensions (timestep, location), " *
            "got $(ndims(dhw))D"
        ))
    end

    if size(dhw, 1) == 0 || size(dhw, 2) == 0
        throw(ArgumentError(
            "DHW YAXArray dimensions must be non-zero, got size $(size(dhw))"
        ))
    end

    return generate_environment(
        convert(Matrix{Float32}, dhw.data); start_year=start_year
    )
end
