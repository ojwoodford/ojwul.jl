export createpackage
export @SR
using PkgTemplates, UserConfig, StaticArrays

function createpackage(name::String)
    template = PkgTemplates.Template(; user=localstring("Github username"), 
                                       license="MIT",
                                       authors=localstring("full name"),
                                       plugins=[Git(), GitHubActions(), ])
    template(name)
end

macro SR(start, finish)
    return :( SVector{$finish-$start+1, Int}($start:$finish) )
end