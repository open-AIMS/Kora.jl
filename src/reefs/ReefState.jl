using Distributions, KernelDensity, Random, StatsBase

using CurveFit

using YAXArrays


struct ReefState{F<:AbstractFloat,   # For Float32
                 P<:Polynomial,      # For Polynomial
                 Y3<:YAXArray{F,3},  # For 3D YAXArrays
                 Y4a<:YAXArray{F,4},
                 Y4b<:YAXArray{F,4},
                 Y4c<:YAXArray{F,4}}  # For 4D YAXArrays
    pop_sample::Y4a
    growth_models::Vector{P}
    survival_models::Vector{P}
    location_scalers::Y3
    dhw_tolerances::Y4b
    mortalities::Y4c
    carrying_capacity::Vector{F}
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

@inline function population_sample(reef_state::ReefState, ts::Int64, loc::Int64)::SubArray{Float32}
    return @inbounds @views reef_state.pop_sample.data[ts, loc, :]
end

function update_sample!(reef_state::ReefState, ts::Int64, loc::Int64, group::Int64, size_dist::AbstractVector{Float32})::Nothing
    @inbounds reef_state.pop_sample.data[ts, loc, group, :] .= size_dist

    return nothing
end

function update_dhw_tol_mean!(reef_state, ts::Int64, grp::Int64, vals::AbstractVector{Float32})::Nothing
    reef_state.dhw_tolerances.data[ts, :, grp, 1] = vals

    return nothing
end

function initialize_reef(; n_timesteps=75, n_locs=100, n_groups=5, sample_size=1000)
    group_names = Symbol.(["group$(i)" for i in 1:n_groups])

    loc_ax = (
        Dim{:scaler}([:growth, :mortality]),
        Dim{:location}(1:n_locs),
        Dim{:group}(group_names)
    )
    location_scalers = YAXArray(loc_ax, ones(Float32, 2, n_locs, n_groups))

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

    # _pop_cache::Matrix{F}        # Temporary store of population
    # _growth_cache::Matrix{F}     # Temporary store of population growth
    # _location_cache::Vector{F}   # Temporary store of population for a given location
    # _pop_buffer::Matrix{F}       # Continually updated store of cover for all locations
    # _location_buffer::Vector{F}  # Continually updated store of cover for all locations
    # _recruit_buffer::Vector{F}   # Continually refreshed store of coral recruits
    reef_state = ReefState(
        size_dist,
        CoralFlow.growth_models,
        CoralFlow.survival_models,
        location_scalers,
        dhw_tol,
        mort,
        fill(1500.0f0, 100),
        zeros(Float32, n_locs, sample_size),
        zeros(Float32, n_locs, sample_size),
        zeros(Float32, n_locs),
        zeros(Float32, n_locs, sample_size*2),
        zeros(Float32, n_locs),
        zeros(Float32, n_locs)
    )

    return reef_state
end

function initialize_coral_population!(
    reef_state::ReefState,
    μ::Float32,
    σ::Float32,
    loc::Int64,
    target_pop_size::Int64;
    max_size::Float32=150.0f0,
    rng::AbstractRNG=Random.GLOBAL_RNG
)::Nothing
    dist = truncated(LogNormal(μ, σ), 1.0f0, max_size)
    sample_sizes = convert.(Float32, rand(rng, dist, n_groups(reef_state), target_pop_size))

    reef_state.dhw_tolerances[1, :, :, At(:mean)] .= 4.0f0
    reef_state.dhw_tolerances[:, :, :, At(:stdev)] .= 2.75f0
    reef_state.pop_sample.data[1, loc, :, :] .= sample_sizes

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
