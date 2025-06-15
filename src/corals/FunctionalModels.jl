abstract type AbstractCoralBehavior end


Base.getindex(x::AbstractCoralBehavior, i::Union{Int,UnitRange}) = x.models[i]
