"""
    bin_edges()::Matrix{Float32}

Helper function defining coral colony diameter bin edges. The values are converted from `cm`
to the desired unit. The default unit is `m`.
"""
function bin_edges()::Matrix{Float32}
    return Matrix{Float32}(
        [
            2.5 5.0 7.5 10.0 20.0 40.0 100.0 150.0;
            2.5 5.0 7.5 10.0 20.0 35.0 50.0 100.0;
            2.5 5.0 7.5 10.0 15.0 20.0 40.0 50.0;
            2.5 5.0 7.5 10.0 20.0 40.0 50.0 100.0;
            2.5 5.0 7.5 10.0 20.0 40.0 50.0 100.0
        ]
    )
end

"""
    diameter_size_classes()::Vector{Matrix{Float32}}

Determine diameter widths for each size class.

See also:
- `bin_edges()`
- `bin_widths()`
"""
function diameter_size_classes()::Vector{Matrix{Float32}}
    edges = bin_edges()
    n_groups = size(edges, 1)
    edges = hcat(zeros(Float32, n_groups), edges)

    # TODO: Use YAXArrays (groups ⋅ class start ⋅ class end)
    return map(x -> [x[1:end-1] x[2:end]], eachrow(edges))
end

"""
    bin_widths()

Helper function defining coral colony diameter bin widths.
"""
function bin_widths()
    return bin_edges()[:, 2:end] .- bin_edges()[:, 1:(end - 1)]
end

"""
    linear_extensions()::Matrix{Float32}

Linear extensions.
"""
function linear_extensions()::Matrix{Float32}
    return [
        0.609456f0 1.07184f0 2.55149f0 5.07988f0 9.45091f0 16.8505f0 8.0f0 0.0f0;
        0.768556f0 1.22085f0 1.86447f0 2.82297f0 3.52938f0 3.00422f0 1.5f0 0.0f0;
        0.190455f0 0.343747f0 0.615467f0 0.97477f0 1.70079f0 2.91729f0 1.45f0 0.0f0;
        0.318034f0 0.47385f0 0.683729f0 0.710587f0 0.581085f0 0.581085f0 0.3f0 0.0f0;
        0.122478f0 0.217702f0 0.382098f0 0.718781f0 1.24172f0 2.08546f0 1.04f0 0.0f0
    ]
end

"""
    survival_rates()::Matrix{Float32}

Survival rates.
"""
function survival_rates()::Matrix{Float32}
    return [
        0.6f0 0.76f0 0.805f0 0.76f0 0.85f0 0.86f0 0.86f0 0.86f0;    # Tabular Acropora
        0.6f0 0.76f0 0.77f0 0.875f0 0.83f0 0.90f0 0.90f0 0.90f0;    # Corymbose Acropora
        0.52f0 0.77f0 0.77f0 0.875f0 0.89f0 0.97621179f0 0.97621179f0 0.97621179f0;                # Corymbose non-Acropora
        0.72f0 0.87f0 0.77f0 0.98f0 0.996931548f0 0.996931548f0 0.996931548f0 0.996931548f0;        # Small massives and encrusting
        0.58f0 0.87f0 0.78f0 0.983568572f0 0.984667677f0 0.984667677f0 0.984667677f0 0.984667677f0  # Large massives
    ]
end
