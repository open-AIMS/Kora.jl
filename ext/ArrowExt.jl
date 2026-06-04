module ArrowExt

using Arrow, DataFrames
using Kora

function Kora._read_arrow_file(filepath::String)::DataFrame
    return DataFrame(Arrow.Table(filepath))
end

end
