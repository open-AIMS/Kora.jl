module CoralFlow

using Random
using Statistics
using Bootstrap

using FLoops

include("./stats.jl")
include("./metrics.jl")
include("./reefs/ReefState.jl")
include("./corals/size_classes.jl")
include("./corals/FunctionalModels.jl")
include("./corals/growth_model.jl")
include("./corals/mortality_model.jl")
include("./corals/recruitment.jl")
include("./reefs/flow_model.jl")
include("./viz/viz.jl")

export
    truncated_standard_normal_mean,
    truncated_normal_mean,
    truncated_normal_cdf

export
    growth_models,
    survival_models
    GrowthModel,
    SurvivalModel

export
    cover_cm_to_m2,
    cover_cm_to_m2!,
    coral_cover,
    recruit_cover

export
    ReefState,
    initialize_reef,
    initialize_coral_population!,
    n_timesteps,
    n_locations,
    n_groups,
    pop_sample_size

export
    population_sample,
    update_sample!

end  # module CoralFlow
