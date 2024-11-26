module CoralFlow

using Random
using Statistics

using FLoops

include("./stats.jl")
include("./metrics.jl")
include("./corals/size_classes.jl")
include("./corals/growth_model.jl")
include("./corals/mort_model.jl")
include("./reefs/ReefState.jl")

export
    truncated_standard_normal_mean,
    truncated_normal_mean,
    truncated_normal_cdf

export
    growth_models,
    survival_models

export
    cover_cm_to_m2,
    cover_cm_to_m2!,
    coral_cover,
    recruit_cover

export
    ReefState,
    initialize_reef,
    initialize_coral_population!

export
    population_sample,
    update_sample!

end  # module CoralFlow
