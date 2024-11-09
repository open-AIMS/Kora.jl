"""
    bin_edges()

Helper function defining coral colony diameter bin edges. The values are converted from `cm`
to the desired unit. The default unit is `m`.
"""
function bin_edges()
    return Matrix(
        [
            2.5 5.0 7.5 10.0 20.0 40.0 100.0 150.0;
            2.5 5.0 7.5 10.0 20.0 35.0 50.0 100.0;
            2.5 5.0 7.5 10.0 15.0 20.0 40.0 50.0;
            2.5 5.0 7.5 10.0 20.0 40.0 50.0 100.0;
            2.5 5.0 7.5 10.0 20.0 40.0 50.0 100.0
        ]
    )
end

function diameter_size_classes()
    edges = bin_edges()
    n_groups = size(edges, 1)
    edges = hcat(zeros(n_groups), edges)

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
    linear_extensions()

Linear extensions.
"""
function linear_extensions()::Matrix{Float64}
    return [
        0.609456 1.07184 2.55149 5.07988 9.45091 16.8505 0.0;
        0.768556 1.22085 1.86447 2.82297 3.52938 3.00422 0.0;
        0.190455 0.343747 0.615467 0.97477 1.70079 2.91729 0.0;
        0.318034 0.47385 0.683729 0.710587 0.581085 0.581085 0.0;
        0.122478 0.217702 0.382098 0.718781 1.24172 2.08546 0.0
    ]
end
