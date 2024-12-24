abstract type AbstractCoral end


Base.getindex(x::AbstractCoral, i::Union{Int, UnitRange}) = x.models[i]

function Base.show(io::IO, ::MIME"text/plain", x::AbstractCoral)
    println(io, "")
    println(io, "RMSE: $(x.rmse_scores)")
    println(io, "R2: $(x.r2_scores)")
end