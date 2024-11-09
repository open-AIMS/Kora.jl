module CoralFlow

using Statistics

include("./stats.jl")
include("./metrics.jl")
include("./corals/size_classes.jl")


export
    truncated_standard_normal_mean,
    truncated_normal_mean,
    truncated_normal_cdf


end  # module CoralFlow
