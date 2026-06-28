using JSON
using Polynomials: Polynomial, coeffs
using Dates: now, DateTime, format, value

const SUPPORTED_FORMAT_VERSIONS = (1,)

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
    tr = obj["train"]
    te = obj["test"]
    train_nt = (
        RMSE     = Float32.(tr["RMSE"]),
        R2       = Float32.(tr["R2"]),
        pearson  = Float32.(tr["pearson"]),
        spearman = Float32.(tr["spearman"]),
        kendall  = Float32.(tr["kendall"])
    )
    test_nt = (
        RMSE     = Float32.(te["RMSE"]),
        R2       = Float32.(te["R2"]),
        pearson  = Float32.(te["pearson"]),
        spearman = Float32.(te["spearman"]),
        kendall  = Float32.(te["kendall"])
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
`"growth"` or `"survival"`, and a `"type"` tag for each model entry of either
`"PolyGrowthFunction"` or `"PolySurvivalFunction"`. An informative error is
raised if any required field is missing or if an unknown type tag is encountered.

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
    growth_fns = PolyGrowthFunction[]
    surv_fns = PolySurvivalFunction[]

    for (i, entry) in enumerate(raw["models"])
        validate_spec(entry, ["type", "name"], "$filepath models[$i]")
        push!(names, String(entry["name"]))
        tag = String(entry["type"])
        if tag == "PolyGrowthFunction"
            push!(growth_fns, _deserialize_poly_growth(entry))
        elseif tag == "PolySurvivalFunction"
            push!(surv_fns, _deserialize_poly_survival(entry))
        else
            error(
                "Unknown model type tag \"$tag\" in \"$filepath\". " *
                "Expected \"PolyGrowthFunction\" or \"PolySurvivalFunction\"."
            )
        end
    end

    performance = _dict_to_performance(raw["performance"])

    if kind == "growth"
        return PolyGrowthModel(names, growth_fns, performance)
    elseif kind == "survival"
        return PolySurvivalModel(names, surv_fns, performance)
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

# Parse "yyyy-mm-ddTHH:MM:SS" by splitting into integer components to avoid
# Dates.parse, which triggers unresolvable dispatch in JuliaC's verifier.
function _parse_iso_datetime(s::String)::DateTime
    return DateTime(
        parse(Int, s[1:4]),   # year
        parse(Int, s[6:7]),   # month
        parse(Int, s[9:10]),  # day
        parse(Int, s[12:13]), # hour
        parse(Int, s[15:16]), # minute
        parse(Int, s[18:19]), # second
    )
end

function check_model_pair_skew(
    growth_path::String, surv_path::String; threshold_seconds::Int=86400
)::Nothing
    function _read_fitted_at(path::String)::Union{String,Nothing}
        content = read(path, String)
        m = match(r"\"fitted_at\"\s*:\s*\"([^\"]+)\"", content)
        m === nothing && return nothing
        cap = m.captures[1]
        cap === nothing && return nothing
        return String(cap::SubString{String})
    end

    ga = _read_fitted_at(growth_path)
    sa = _read_fitted_at(surv_path)

    if isnothing(ga) || isnothing(sa)
        @warn "check_model_pair_skew: one or both files missing fitted_at." growth =
            growth_path survival = surv_path
        return nothing
    end

    try
        gt = _parse_iso_datetime(ga)
        st = _parse_iso_datetime(sa)
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
