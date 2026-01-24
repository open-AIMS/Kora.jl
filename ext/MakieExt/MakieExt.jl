module MakieExt

using Random
using Distributions
using OrderedCollections

using YAXArrays
using Makie
using Makie: distinguishable_colors
using CoralFlow
using CoralFlow.Bootstrap

const SUPPORTED_FILETYPES = ("png", "svg", "gif", "jpg", "jpeg")

COLOR_MAP = :Paired_12
FGROUP_COLOR = distinguishable_colors(8)[3:end]
FLABELS = [
    "Tabular Acropora", "Corymbose Acropora",
    "branching non-Acropora", "Small massives", "Large massives"
]

function _display_or_save(
    fig::Figure, type_name::String, group_name::String, save_path::Union{String,Nothing}
)::Nothing
    # Display or save
    if !isnothing(save_path)
        if any(endswith.(Ref(save_path), SUPPORTED_FILETYPES))
            save(save_path, fig)
        elseif isdir(save_path)
            save(joinpath(save_path, "$(type_name)_$(group_name).png"), fig)
        else
            msg = """
            Provided figure savepath is not an existing directory or a supported filetype.
            If a file extension is specified, it must be one of: $(SUPPORTED_FILETYPES)
            """
            throw(ArgumentError(msg))
        end

        return nothing
    end

    display(fig)

    return nothing
end

include("performance.jl")
include("run_analysis.jl")
include("ensemble_analysis.jl")

end