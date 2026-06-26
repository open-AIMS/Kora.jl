const TARGET_GROUPS = CSV.read(
    joinpath(_kora_assets_dir(), "target_groups.csv"), DataFrame; types=String
)[
    :, :functional_group
]

const GROUP_NAMES = CSV.read(
    joinpath(_kora_assets_dir(), "target_groups.csv"), DataFrame; types=String
)[
    :, :group_name
]

include("FunctionalModels.jl")
include("growth_model.jl")
include("mortality_model.jl")
include("recruitment.jl")
include("size_classes.jl")
