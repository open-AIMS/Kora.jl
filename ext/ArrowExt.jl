module ArrowExt

using Arrow, DataFrames
using Kora

function __init__()
    Kora._FILE_READERS[".arrow"] = (f) -> DataFrame(Arrow.Table(f))
end

end
