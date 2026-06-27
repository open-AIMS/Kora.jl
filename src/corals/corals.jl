const TARGET_GROUPS = [
    "acro_table",
    "acro_corym",
    "corym_non_acro",
    "small_massive",
    "large_massive"
]

const GROUP_NAMES = [
    "tabular Acropora",
    "corymbose Acropora",
    "branching non-Acropora",
    "Small Massives",
    "Large Massives"
]

include("FunctionalModels.jl")
include("growth_model.jl")
include("mortality_model.jl")
include("recruitment.jl")
include("size_classes.jl")
