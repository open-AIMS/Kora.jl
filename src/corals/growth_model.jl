using CurveFit

# Private consts - only used in this file
const _euler_f32b = Float32(ℯ)


"""
    growth_inflection_point()::Vector{Float32}

Value indicates when to begin adjusting growth rates as cover approaches the maximum
available hard substrate.

See: `space_constraint()`
"""
function growth_inflection_point()::Vector{Float32}
    return Float32[0.96, 0.96, 0.96, 0.96, 0.96, 0.96]
end

"""
    space_constraint(x, k; x0=0.96)

Modify growth as coral cover approaches total habitable area.

# Arguments
- `x` : proportion of available space currently taken
- `k` : steepness parameter, where the larger the value, the faster the drop (suggest k=20).
- `x0` : Inflection point where growth begins to be constrained (default: 96% of available area).
"""
@inline function space_constraint(x::F, k::F; x0::F=0.96f0)::F where {F<:Float32}
    if x >= 1.0
        return 0.0
    end

    return 1.0f0 / (1.0f0 + _euler_f32b^(k*(x - x0)))
end

# Growth models need to account for available space... (see space constraint)
@inline function growth(model::M, diam::F, reef_cover::F, grp_mod::F, loc_scaler::F)::F where {M,F<:Float32}
    return diam + (model(diam) - diam) * space_constraint(reef_cover, 20.0f0; x0=grp_mod) * loc_scaler
end
function growth!(model::M, diam::AbstractMatrix{F}, reef_cover::Vector{F}, grp_mod::F, loc_scalers::Vector{F})::Nothing where {M,F<:Float32}
    Threads.@threads for i in eachindex(reef_cover)
        # Combine model evaluation, max, and scaling in one pass
        @inbounds @views diam[i, :][diam[i, :] .!= 0.0f0] = growth.(model, diam[i, :][diam[i, :] .!= 0.0f0], reef_cover[i], grp_mod, loc_scalers[i])
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
    _pop = pop[loc, :][pop[loc, :] .> 0.0]
    next_pop = @view(ecostate._pop_cache[loc, :])
    max_sample_size = size(ecostate.pop_sample, 4)

    # If population is below max sample size, keep everything
    if length(_pop) <= max_sample_size
        next_pop[1:length(_pop)] .= _pop
        if length(_pop) < max_sample_size
            next_pop[(length(_pop)+1):end] .= 0.0f0
        end
    else
        # Only resample if we exceed max sample size
        # Use size-stratified sampling to maintain size distribution
        size_weights = zeros(Float32, length(_pop))
        for (sz_end, sz_start) in eachrow(class_diams)
            in_class = sz_start .< _pop .<= sz_end
            if any(in_class)
                # Weight by colony size to maintain proper cover representation
                size_weights[in_class] .= _pop[in_class] ./ sum(_pop[in_class])
            end
        end

        # Sample proportionally to maintain size structure
        sampled_idx = sample(rng, 1:length(_pop), Weights(size_weights), max_sample_size, replace=false)
        next_pop .= _pop[sampled_idx]
    end

    update_sample!(ecostate, ts, loc, grp, next_pop)

    return nothing
end

growth_models = Function[]
growth_model_coefs = Vector{Float32}[]
growth_rmse_scores = Float32[]
growth_r2_scores = Float32[]
for (xi, yi) in zip(eachrow(CoralFlow.bin_edges()[:, 1:end-1]), eachrow(CoralFlow.linear_extensions()))
    m = curve_fit(Polynomial, xi, yi, 3)
    model = x -> x .< xi[1] ? yi[1] : max(m(x), 0.1)

    push!(growth_rmse_scores, CoralFlow.RMSE(model.(xi), yi))
    push!(growth_r2_scores, CoralFlow.R2(model.(xi), vec(yi)))
    push!(growth_models, model)

    # Annoyingly, different regression types have different fieldnames.
    # Although we've settled on `Polynomial` for now, leaving this try/catch in here
    # in case this changes.
    try
        push!(growth_model_coefs, m.coefs)
    catch err
        push!(growth_model_coefs, m.coeffs)
    end
end

@kwdef struct GrowthModel <: AbstractCoral
    models::Vector{Function}
    rmse_scores::Vector{Float64} = []
    r2_scores::Vector{Float64} = []
end

function GrowthModel(models)
    return GrowthModel(models, Float32[], Float32[])
end
