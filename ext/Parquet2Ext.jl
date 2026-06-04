module Parquet2Ext

using Parquet2, DataFrames
using Kora

function __init__()
    Kora._FILE_READERS[".parquet"] = (f) -> DataFrame(Parquet2.readfile(f))
end

end
