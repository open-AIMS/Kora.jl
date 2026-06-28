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
    content = read(filepath, String)

    # --- format_version ---
    fv_m = match(r"\"format_version\"\s*:\s*(\d+)", content)
    fv_m === nothing && error("missing format_version in $filepath")
    fvc = fv_m.captures[1]
    fvc === nothing && error("format_version capture failed in $filepath")
    fv = parse(Int, fvc::SubString{String})
    fv == 1 || error("unsupported format_version $fv in $filepath")

    # --- model_kind ---
    mk_m = match(r"\"model_kind\"\s*:\s*\"([^\"]+)\"", content)
    mk_m === nothing && error("missing model_kind in $filepath")
    mkc = mk_m.captures[1]
    mkc === nothing && error("model_kind capture failed in $filepath")
    model_kind = String(mkc::SubString{String})

    # --- models array: extract each {...} block by balanced-brace scan ---
    models_start = findfirst("\"models\"", content)
    models_start === nothing && error("missing models array in $filepath")
    arr_start = findnext('[', content, last(models_start))
    arr_start === nothing && error("missing models [ in $filepath")

    names = String[]
    growth_fns   = PolyGrowthFunction[]
    survival_fns = PolySurvivalFunction[]

    pos = arr_start + 1
    while pos <= length(content)
        obj_start = findnext('{', content, pos)
        obj_start === nothing && break
        depth = 0
        obj_end = obj_start
        for i in obj_start:length(content)
            c = content[i]
            if c == '{'
                depth += 1
            elseif c == '}'
                depth -= 1
                if depth == 0
                    obj_end = i
                    break
                end
            end
        end
        block = content[obj_start:obj_end]

        # Extract dtype
        dt_m = match(r"\"dtype\"\s*:\s*\"([^\"]+)\"", block)
        dt_m === nothing && error("missing dtype in model block")
        dtc = dt_m.captures[1]
        dtc === nothing && error("dtype capture failed")
        dtype_str = String(dtc::SubString{String})

        # Extract scalar fields
        function extract_str(pat, blk)
            mm = match(pat, blk)
            mm === nothing && error("pattern $pat not found in block")
            cc = mm.captures[1]
            cc === nothing && error("capture failed for $pat")
            return String(cc::SubString{String})
        end

        entry_name  = extract_str(r"\"name\"\s*:\s*\"([^\"]+)\"",  block)
        entry_type  = extract_str(r"\"type\"\s*:\s*\"([^\"]+)\"",  block)

        # Extract numeric scalars as strings for parsing
        mn_x_m = match(r"\"min_x\"\s*:\s*([-\d.eE+]+)", block)
        mn_y_m = match(r"\"min_y\"\s*:\s*([-\d.eE+]+)", block)
        mx_x_m = match(r"\"max_x\"\s*:\s*([-\d.eE+]+)", block)
        mx_y_m = match(r"\"max_y\"\s*:\s*([-\d.eE+]+)", block)
        (mn_x_m === nothing || mn_y_m === nothing || mx_x_m === nothing || mx_y_m === nothing) &&
            error("missing min/max fields in block")
        mnxc = mn_x_m.captures[1]; mnxc === nothing && error("min_x capture")
        mnyc = mn_y_m.captures[1]; mnyc === nothing && error("min_y capture")
        mxxc = mx_x_m.captures[1]; mxxc === nothing && error("max_x capture")
        mxyc = mx_y_m.captures[1]; mxyc === nothing && error("max_y capture")
        min_x_str = mnxc::SubString{String}
        min_y_str = mnyc::SubString{String}
        max_x_str = mxxc::SubString{String}
        max_y_str = mxyc::SubString{String}

        # Extract poly_coeffs array content
        pc_m = match(r"\"poly_coeffs\"\s*:\s*\[([^\]]*)\]", block)
        pc_m === nothing && error("missing poly_coeffs in block")
        pcc = pc_m.captures[1]
        pcc === nothing && error("poly_coeffs capture failed")
        coeffs_substr = pcc::SubString{String}

        push!(names, entry_name)

        if dtype_str == "float32"
            min_x = Float32(parse(Float64, min_x_str))
            min_y = Float32(parse(Float64, min_y_str))
            max_x = Float32(parse(Float64, max_x_str))
            max_y = Float32(parse(Float64, max_y_str))
            coeff_matches = eachmatch(r"[-\d.eE+]+", coeffs_substr)
            coeffs = Float32[Float32(parse(Float64, cm.match)) for cm in coeff_matches]
            if entry_type == "PolyGrowthFunction" || model_kind == "growth"
                push!(growth_fns, PolyGrowthFunction(min_x, min_y, max_x, max_y, Polynomial(coeffs)))
            elseif entry_type == "PolySurvivalFunction" || model_kind == "survival"
                push!(survival_fns, PolySurvivalFunction(min_x, min_y, max_x, max_y, Polynomial(coeffs)))
            else
                error("unknown model type $entry_type")
            end
        elseif dtype_str == "float64"
            min_x = parse(Float64, min_x_str)
            min_y = parse(Float64, min_y_str)
            max_x = parse(Float64, max_x_str)
            max_y = parse(Float64, max_y_str)
            coeff_matches = eachmatch(r"[-\d.eE+]+", coeffs_substr)
            coeffs = Float64[parse(Float64, cm.match) for cm in coeff_matches]
            if entry_type == "PolyGrowthFunction" || model_kind == "growth"
                push!(growth_fns, PolyGrowthFunction(min_x, min_y, max_x, max_y, Polynomial(coeffs)))
            elseif entry_type == "PolySurvivalFunction" || model_kind == "survival"
                push!(survival_fns, PolySurvivalFunction(min_x, min_y, max_x, max_y, Polynomial(coeffs)))
            else
                error("unknown model type $entry_type")
            end
        else
            error("unknown dtype $dtype_str in $filepath")
        end

        pos = obj_end + 1
    end

    # --- performance section ---
    perf_start = findfirst("\"performance\"", content)
    perf_start === nothing && error("missing performance section in $filepath")

    function extract_metric(metric_name, split_name)
        # Find the split section (train/test) then the metric array
        split_pat = Regex("\"$split_name\"\\s*:\\s*\\{")
        sm_m = match(split_pat, content)
        sm_m === nothing && error("missing $split_name in performance")
        # Find metric array within the region after perf_start
        search_region = content[first(perf_start):end]
        split_m2 = match(split_pat, search_region)
        split_m2 === nothing && error("missing $split_name in performance region")
        split_pos = first(perf_start) + split_m2.offset - 1
        met_pat = Regex("\"$metric_name\"\\s*:\\s*\\[([^\\]]*?)\\]")
        # search from split_pos
        met_m = match(met_pat, content[split_pos:end])
        met_m === nothing && error("missing $metric_name in $split_name")
        arr_c = met_m.captures[1]
        arr_c === nothing && error("$metric_name capture failed")
        arr_str = arr_c::SubString{String}
        vals = Float32[]
        for vm in eachmatch(r"[-\d.eE+]+|null", arr_str)
            push!(vals, vm.match == "null" ? NaN32 : Float32(parse(Float64, vm.match)))
        end
        return vals
    end

    train_rmse    = extract_metric("RMSE",     "train")
    train_r2      = extract_metric("R2",       "train")
    train_pearson = extract_metric("pearson",  "train")
    train_spearman= extract_metric("spearman", "train")
    train_kendall = extract_metric("kendall",  "train")

    test_rmse     = extract_metric("RMSE",     "test")
    test_r2       = extract_metric("R2",       "test")
    test_pearson  = extract_metric("pearson",  "test")
    test_spearman = extract_metric("spearman", "test")
    test_kendall  = extract_metric("kendall",  "test")

    train_perf = (RMSE=train_rmse, R2=train_r2, pearson=train_pearson,
                  spearman=train_spearman, kendall=train_kendall)
    test_perf  = (RMSE=test_rmse,  R2=test_r2,  pearson=test_pearson,
                  spearman=test_spearman,  kendall=test_kendall)
    performance = (train=train_perf, test=test_perf)

    if model_kind == "growth"
        isempty(growth_fns) && error("no growth functions parsed from $filepath")
        return PolyGrowthModel(names, growth_fns, performance)
    elseif model_kind == "survival"
        isempty(survival_fns) && error("no survival functions parsed from $filepath")
        return PolySurvivalModel(names, survival_fns, performance)
    else
        error("unknown model_kind $model_kind in $filepath")
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
