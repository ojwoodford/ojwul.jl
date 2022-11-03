export createpackage, SR, valuedispatch1to32
using PkgTemplates, UserConfig, StaticArrays

function valuedispatch1to32(fun, val)
    if val <= 16
        if val <= 8
            if val <= 4
                if val <= 2
                    if val == 2
                        return fun(Val(2))
                    else
                        return fun(Val(1))
                    end
                else
                    if val == 4
                        return fun(Val(4))
                    else
                        return fun(Val(3))
                    end
                end
            else
                if val <= 6
                    if val == 6
                        return fun(Val(6))
                    else
                        return fun(Val(5))
                    end
                else
                    if val == 8
                        return fun(Val(8))
                    else
                        return fun(Val(7))
                    end
                end
            end
        else
            if val <= 12
                if val <= 10
                    if val == 10
                        return fun(Val(10))
                    else
                        return fun(Val(11))
                    end
                else
                    if val == 12
                        return fun(Val(12))
                    else
                        return fun(Val(11))
                    end
                end
            else
                if val <= 14
                    if val == 14
                        return fun(Val(14))
                    else
                        return fun(Val(13))
                    end
                else
                    if val == 16
                        return fun(Val(16))
                    else
                        return fun(Val(15))
                    end
                end
            end
        end
    else
        if val <= 24
            if val <= 20
                if val <= 18
                    if val == 18
                        return fun(Val(18))
                    else
                        return fun(Val(17))
                    end
                else
                    if val == 20
                        return fun(Val(20))
                    else
                        return fun(Val(19))
                    end
                end
            else
                if val <= 22
                    if val == 22
                        return fun(Val(22))
                    else
                        return fun(Val(21))
                    end
                else
                    if val == 24
                        return fun(Val(24))
                    else
                        return fun(Val(23))
                    end
                end
            end
        else
            if val <= 28
                if val <= 26
                    if val == 26
                        return fun(Val(26))
                    else
                        return fun(Val(25))
                    end
                else
                    if val == 28
                        return fun(Val(28))
                    else
                        return fun(Val(27))
                    end
                end
            else
                if val <= 30
                    if val == 30
                        return fun(Val(30))
                    else
                        return fun(Val(29))
                    end
                else
                    if val == 32
                        return fun(Val(32))
                    else
                        return fun(Val(31))
                    end
                end
            end
        end
    end
end

function createpackage(name::String)
    template = PkgTemplates.Template(; user=localstring("Github username"), 
                                       license="MIT",
                                       authors=localstring("full name"),
                                       plugins=[Git(), GitHubActions(), ])
    template(name)
end

const SR = StaticArrays.SUnitRange