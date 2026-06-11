using JSON
using Polynomials: Polynomial, coeffs
using Dates: now, DateTime, format, value

const SUPPORTED_FORMAT_VERSIONS = (1,)

# ---------------------------------------------------------------------------
# Registry
# Deserializer functions are defined at include-time (stable named functions,
# no closures). register_model_type! is called from __init__ so no function
# reference is baked into the precompile cache.
#
# Thread safety: _MODEL_TYPE_REGISTRY is populated exclusively in __init__ and
# is read-only after startup. Concurrent registration from user-spawned tasks
# is not supported and is not thread-safe.
# ---------------------------------------------------------------------------

const _MODEL_TYPE_REGISTRY = Dict{String,Function}()

"""
    register_model_type!(tag::String, deserializer::Function)::Nothing

Register a custom model type so that `load_models` can deserialise it from JSON.

The `tag` string must match the `"type"` field written by the custom serialiser into
each model entry. The `deserializer` receives the raw parsed JSON `Dict` for
that entry and must return a callable representing the fitted model function.

The built-in types `"PolyGrowthFunction"` and `"PolySurvivalFunction"` are
registered automatically when the package is loaded via `__init__`. Registration
is not thread-safe; all calls must occur before concurrent model loading begins.

# Arguments
- `tag::String` : Unique string identifier stored in the JSON file.
- `deserializer::Function` : Function with signature
  `(entry::AbstractDict) -> Function` that reconstructs a model callable from
  the serialised `Dict`.

# Returns
`Nothing`

# See Also
[`load_models`](@ref), [`save_models`](@ref)
"""
function register_model_type!(tag::String, deserializer::Function)::Nothing
    _MODEL_TYPE_REGISTRY[tag] = deserializer
    return nothing
end

# Helpers

function validate_spec(spec::AbstractDict, required::Vector{String}, path::String)::Nothing
    missing_fields = filter(f -> !haskey(spec, f), required)
    if !isempty(missing_fields)
        error(
            "JSON at \"$path\" is missing required field(s): " * join(missing_fields, ", ")
        )
    end
    return nothing
end

function _checked_type_lookup(tag::String)::Function
    if !haskey(_MODEL_TYPE_REGISTRY, tag)
        valid = join(sort(collect(keys(_MODEL_TYPE_REGISTRY))), ", ")
        error("Unknown model type tag \"$tag\". Valid tags: $valid")
    end
    return _MODEL_TYPE_REGISTRY[tag]
end

_dtype_str(::Type{Float32}) = "float32"
_dtype_str(::Type{Float64}) = "float64"
_dtype_str(::Type{T}) where {T} = lowercase(string(T))

function _dtype_type(s::AbstractString)
    s == "float32" && return Float32
    s == "float64" && return Float64
    return error("Unknown dtype \"$s\". Expected \"float32\" or \"float64\".")
end

const _METRIC_NAMES = String.(Symbol.(ALL_METRICS))

function _performance_to_dict(perf::NamedTuple)::Dict{String,Any}
    _sanitize(v) = isfinite(v) ? v : nothing
    return Dict{String,Any}(
        "train" => Dict{String,Any}(
            k => _sanitize.(collect(getfield(perf.train, Symbol(k)))) for k in _METRIC_NAMES
        ),
        "test" => Dict{String,Any}(
            k => _sanitize.(collect(getfield(perf.test, Symbol(k)))) for k in _METRIC_NAMES
        )
    )
end

function _dict_to_performance(obj::AbstractDict)::NamedTuple
    train_nt = NamedTuple(
        Symbol(k) => Float32.(obj["train"][k]) for k in _METRIC_NAMES
    )
    test_nt = NamedTuple(
        Symbol(k) => Float32.(obj["test"][k]) for k in _METRIC_NAMES
    )
    return (train=train_nt, test=test_nt)
end

# Deserializer functions (pure input->output, no module-global captures)

function _deserialize_poly_growth(entry::AbstractDict)::PolyGrowthFunction
    validate_spec(
        entry,
        ["dtype", "min_x", "min_y", "max_x", "max_y", "poly_coeffs"],
        "models[type=PolyGrowthFunction]"
    )
    T = _dtype_type(entry["dtype"])
    return PolyGrowthFunction(
        T(entry["min_x"]),
        T(entry["min_y"]),
        T(entry["max_x"]),
        T(entry["max_y"]),
        Polynomial(T.(entry["poly_coeffs"]))
    )
end

function _deserialize_poly_survival(entry::AbstractDict)::PolySurvivalFunction
    validate_spec(
        entry,
        ["dtype", "min_x", "min_y", "max_x", "max_y", "poly_coeffs"],
        "models[type=PolySurvivalFunction]"
    )
    T = _dtype_type(entry["dtype"])
    return PolySurvivalFunction(
        T(entry["min_x"]),
        T(entry["min_y"]),
        T(entry["max_x"]),
        T(entry["max_y"]),
        Polynomial(T.(entry["poly_coeffs"]))
    )
end

# Serialization helpers

function _model_to_entry(f::PolyGrowthFunction{T}, name::String)::Dict{String,Any} where {T}
    return Dict{String,Any}(
        "type" => "PolyGrowthFunction",
        "name" => name,
        "dtype" => _dtype_str(T),
        "min_x" => f.min_x,
        "min_y" => f.min_y,
        "max_x" => f.max_x,
        "max_y" => f.max_y,
        "poly_coeffs" => collect(coeffs(f.poly))
    )
end

function _model_to_entry(
    f::PolySurvivalFunction{T}, name::String
)::Dict{String,Any} where {T}
    return Dict{String,Any}(
        "type" => "PolySurvivalFunction",
        "name" => name,
        "dtype" => _dtype_str(T),
        "min_x" => f.min_x,
        "min_y" => f.min_y,
        "max_x" => f.max_x,
        "max_y" => f.max_y,
        "poly_coeffs" => collect(coeffs(f.poly))
    )
end

# ---------------------------------------------------------------------------
# save_models
# ---------------------------------------------------------------------------

"""
    save_models(m::PolyGrowthModel, filepath::String; region::String="")::Nothing
    save_models(m::PolySurvivalModel, filepath::String; region::String="")::Nothing

Serialise a fitted model collection to a versioned JSON file at `filepath`.

The file records model kind (`"growth"` or `"survival"`), fit timestamp,
polynomial coefficients, numeric range limits, and per-metric train/test
performance. It can be reloaded without information loss using `load_models`.

# Arguments
- `m` : Fitted model collection to serialise (`PolyGrowthModel` or
  `PolySurvivalModel`).
- `filepath::String` : Destination path including filename and `.json` extension.
- `region::String` : Optional label stored in the file header for reference.
  Does not affect the model data (default: `""`).

# Returns
`Nothing`

# See Also
[`load_models`](@ref), [`register_model_type!`](@ref),
[`fit_growth_models`](@ref), [`fit_survival_models`](@ref)
"""
function save_models(
    m::PolyGrowthModel, filepath::String; region::String=""
)::Nothing
    entries = [_model_to_entry(m.models[i], m.names[i]) for i in eachindex(m.names)]
    doc = Dict{String,Any}(
        "format_version" => 1,
        "model_kind" => "growth",
        "region" => region,
        "fitted_at" => format(now(), "yyyy-mm-ddTHH:MM:SS"),
        "models" => entries,
        "performance" => _performance_to_dict(m.performance)
    )
    open(filepath, "w") do io
        JSON.print(io, doc, 2)
    end
    return nothing
end

function save_models(
    m::PolySurvivalModel, filepath::String; region::String=""
)::Nothing
    entries = [_model_to_entry(m.models[i], m.names[i]) for i in eachindex(m.names)]
    doc = Dict{String,Any}(
        "format_version" => 1,
        "model_kind" => "survival",
        "region" => region,
        "fitted_at" => format(now(), "yyyy-mm-ddTHH:MM:SS"),
        "models" => entries,
        "performance" => _performance_to_dict(m.performance)
    )
    open(filepath, "w") do io
        JSON.print(io, doc, 2)
    end
    return nothing
end

# ---------------------------------------------------------------------------
# load_models
# ---------------------------------------------------------------------------

"""
    load_models(filepath::String)::Union{PolyGrowthModel, PolySurvivalModel}

Deserialise a model collection from a versioned JSON file previously written by
`save_models`.

The file must declare a supported `format_version`, a `model_kind` of either
`"growth"` or `"survival"`, and a `"type"` tag for each model entry that has
been registered via `register_model_type!`. An informative error is raised if
any required field is missing or if an unknown type tag is encountered.

# Arguments
- `filepath::String` : Path to the JSON model file.

# Returns
`Union{PolyGrowthModel, PolySurvivalModel}` : The deserialised model collection,
ready for use as the `growth_models` or `survival_models` argument to
`initialize_reef`.

# Examples
```jldoctest
julia> using Kora

julia> path = joinpath(pkgdir(Kora), "assets", "models",
           "offshore_north_growth_models.json");

julia> m = load_models(path);

julia> typeof(m)
Kora.PolyGrowthModel
```

# See Also
[`save_models`](@ref), [`register_model_type!`](@ref),
[`initialize_reef`](@ref)
"""
function load_models(filepath::String)::Union{PolyGrowthModel,PolySurvivalModel}
    raw = JSON.parse(read(filepath, String))

    validate_spec(
        raw,
        ["format_version", "model_kind", "fitted_at", "models", "performance"],
        filepath
    )

    fv = raw["format_version"]
    if fv ∉ SUPPORTED_FORMAT_VERSIONS
        error(
            "Unsupported format_version=$fv in \"$filepath\". " *
            "Supported: $(SUPPORTED_FORMAT_VERSIONS)"
        )
    end

    kind = raw["model_kind"]
    names = String[]
    models = Function[]

    for (i, entry) in enumerate(raw["models"])
        validate_spec(entry, ["type", "name"], "$filepath models[$i]")
        type_tag = String(entry["type"])
        deserializer = _checked_type_lookup(type_tag)
        push!(names, String(entry["name"]))
        push!(models, deserializer(entry))
    end

    performance = _dict_to_performance(raw["performance"])

    if kind == "growth"
        return PolyGrowthModel(names, models, performance)
    elseif kind == "survival"
        return PolySurvivalModel(names, models, performance)
    else
        error(
            "Unknown model_kind=\"$kind\" in \"$filepath\". " *
            "Expected \"growth\" or \"survival\"."
        )
    end
end

# ---------------------------------------------------------------------------
# Checked accessors (issue 7: fail with context rather than silently at call site)
# ---------------------------------------------------------------------------

function get_growth_models()::PolyGrowthModel
    global growth_models
    if isnothing(growth_models)
        error(
            "Growth models are not loaded. Check for warnings emitted during `using Kora`."
        )
    end
    return growth_models::PolyGrowthModel
end

function get_survival_models()::PolySurvivalModel
    global survival_models
    if isnothing(survival_models)
        error(
            "Survival models are not loaded. Check for warnings emitted during `using Kora`."
        )
    end
    return survival_models::PolySurvivalModel
end

# ---------------------------------------------------------------------------
# Timestamp skew check (issue 11)
# ---------------------------------------------------------------------------

function check_model_pair_skew(
    growth_path::String, surv_path::String; threshold_seconds::Int=86400
)::Nothing
    function _read_fitted_at(path::String)::Union{String,Nothing}
        raw = JSON.parse(read(path, String))
        return haskey(raw, "fitted_at") ? String(raw["fitted_at"]) : nothing
    end

    ga = _read_fitted_at(growth_path)
    sa = _read_fitted_at(surv_path)

    if isnothing(ga) || isnothing(sa)
        @warn "check_model_pair_skew: one or both files missing fitted_at." growth =
            growth_path survival = surv_path
        return nothing
    end

    fmt = "yyyy-mm-ddTHH:MM:SS"
    try
        gt = DateTime(ga, fmt)
        st = DateTime(sa, fmt)
        diff_s = abs(value(gt - st)) / 1_000
        if diff_s > threshold_seconds
            @warn(
                "Growth and survival model timestamps differ by $(round(Int, diff_s))s " *
                    "(threshold: $(threshold_seconds)s). " *
                    "They may have been fitted from different datasets.",
                growth_fitted_at = ga,
                survival_fitted_at = sa,
                growth_path = growth_path,
                survival_path = surv_path
            )
        end
    catch e
        @warn "check_model_pair_skew: could not parse timestamps." growth = ga survival = sa exception =
            e
    end

    return nothing
end
