module Parquet2Ext

using Parquet2, DataFrames
using Kora

function Kora._read_parquet_file(filepath::String)::DataFrame
    return DataFrame(Parquet2.readfile(filepath))
end

end
