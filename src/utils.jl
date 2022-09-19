export createpackage, SR
using PkgTemplates, UserConfig, StaticArrays

function createpackage(name::String)
    template = PkgTemplates.Template(; user=localstring("Github username"), 
                                       license="MIT",
                                       authors=localstring("full name"),
                                       plugins=[Git(), GitHubActions(), ])
    template(name)
end

const SR = StaticArrays.SUnitRange