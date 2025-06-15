const TARGET_GROUPS = CSV.read(
    joinpath(ASSET_DIR, "target_groups.csv"), DataFrame; types=String
)[:, :functional_group]

include("FunctionalModels.jl")
include("growth_model.jl")
include("mortality_model.jl")
include("recruitment.jl")
include("size_classes.jl")
