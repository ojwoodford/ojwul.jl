export Robustifier, NoRobust, HuberKernel, GemanMcclureKernel
export robustify, robustifyarray
import Plots.plot

# Robustification
abstract type Robustifier end

function robustifyarray(kernel::Robustifier, cost::AbstractArray)
    weight = similar(cost)
    robcost = similar(cost)
    for ind in eachindex(cost)
        robcost[ind], weight[ind] = robustify(kernel, cost[ind])
    end
    return robcost, weight
end


struct NoRobust{T<:Real} <: Robustifier
    height::T
    height_sqrt::T
end

NoRobust(h) = HuberKernel(h, sqrt(h))

function robustify(kernel::NoRobust, cost::Number)
    return cost * kernel.height, kernel.height_sqrt
end


struct HuberKernel{T<:Real} <: Robustifier
    width::T
    width_squared::T
    height::T
    height_sqrt::T
end

HuberKernel(w, h) = HuberKernel(w, w*w, h, sqrt(h))

function robustify(kernel::HuberKernel, cost::Number)
    if cost < kernel.width_squared
        return cost * kernel.height, kernel.height_sqrt
    end
    sqrtcost = sqrt(cost)
    return (sqrtcost * (kernel.width * 2) - kernel.width_squared) * kernel.height, kernel.width * kernel.height_sqrt / sqrtcost
end


struct GemanMcclureKernel{T<:Real} <: Robustifier
    width_squared::T
    height::T
    height_sqrt::T
end

GemanMcclureKernel(w, h) = GemanMcclureKernel(w*w, h/(w*w), sqrt(h)/w)

function robustify(kernel::GemanMcclureKernel, cost::Number)
    w = kernel.width_squared / (cost + kernel.width_squared)
    return cost * w * kernel.height, w * w * kernel.height_sqrt
end


function displaykernel(kernel, maxval=1)
    x = LinRange(0, maxval, 1000)
    cost, weight = robustifyarray(kernel, x .^ 2)
    plot(x, [cost, weight])
end
