export filehash
using SHA

function filehash(fpath::String)
open(fpath) do file
    return bytes2hex(sha256(file))
end
end