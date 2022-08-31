export AbstractCamera, SimpleCamera, NoDistortionCamera, ExtendedUnifiedCamera
export ideal2image, image2ideal, pixel2image, image2pixel, update, nvars
import ojwul.AbstractVariable
abstract type AbstractCamera <: AbstractVariable end
using StaticArrays, LinearAlgebra


function image2pixel(halfimsz, x)
    return x .* halfimsz[1] .+ halfimsz
end

function pixel2image(halfimsz, x)
    return (x .- halfimsz) ./ halfimsz[1]
end

function pixel2image(halfimsz, x, W)
    return ((x .- halfimsz) ./ halfimsz[1], W .* halfimsz[1])
end


struct SimpleCamera{T<:Number} <: AbstractCamera
    f::T
end
function nvars(var::SimpleCamera)
    return 1
end
function update(var::SimpleCamera, updatevec)
    return SimpleCamera(var.f * exp(updatevec[1]))
end

function ideal2image(camera::SimpleCamera, x)
    return x .* camera.f
end

function image2ideal(camera::SimpleCamera, x)
    return x ./ camera.f
end

function image2ideal(camera::SimpleCamera, x, W)
    return x ./ camera.f, W .* camera.f
end


struct NoDistortionCamera{T<:Number} <: AbstractCamera
    f::SVector{2, T}
    c::SVector{2, T}
end
function nvars(var::NoDistortionCamera)
    return 4
end
function update(var::NoDistortionCamera, updatevec)
    return NoDistortionCamera(SVector(var.f[1] * exp(updatevec[1]), var.f[2] * exp(updatevec[2])), 
                              SVector(var.c[1] + updatevec[3], var.c[2] + updatevec[4]))
end

function ideal2image(camera::NoDistortionCamera, x)
    return x .* camera.f .+ camera.c
end

function image2ideal(camera::NoDistortionCamera, x)
    return (x .- camera.c) ./ camera.f
end

function image2ideal(camera::NoDistortionCamera, x, W)
    return (x .- camera.c) ./ camera.f, W .* camera.f' 
end


struct ExtendedUnifiedCamera{T<:Number} <: AbstractCamera
    f::SVector{2, T}
    c::SVector{2, T}
    alpha::T
    beta::T
end
function nvars(var::ExtendedUnifiedCamera)
    return 6
end
function update(var::ExtendedUnifiedCamera, updatevec)
    alpha = var.alpha * exp(updatevec[5])
    alpha = alpha / (1 + (alpha - var.alpha))
    beta = (var.beta != 0 ? var.beta : eps(var.beta)) * exp(updatevec[6])
    return ExtendedUnifiedCamera(SVector(var.f[1] * exp(updatevec[1]), var.f[2] * exp(updatevec[2])), 
                                 SVector(var.c[1] + updatevec[3], var.c[2] + updatevec[4]), alpha, beta)
end

function ideal2image(camera::ExtendedUnifiedCamera, x)
    z = 1 ./ (1 .+ camera.alpha .* (sqrt.(camera.beta .* sum(abs2, x, dims=1) .+ 1) .- 1))
    return (x .* z) .* camera.f .+ camera.c
end

function image2ideal(camera::ExtendedUnifiedCamera, x)
    y = (x .- camera.c) ./ camera.f
    z = camera.beta .* sum(abs2, y, dims=1)
    z = (1 .+ camera.alpha .* (sqrt.(1 .+ z .* (1 .- 2 .* camera.alpha)) .- 1)) ./ (1 .- z .* (camera.alpha ^ 2))
    return y .* z
end

function image2ideal(camera::ExtendedUnifiedCamera, x, W)
    y = (x - camera.c) ./ camera.f
    t = 1 - 2 * camera.alpha
    n = y' * y
    u = camera.beta * n
    v = 1 - u * (camera.alpha ^ 2)
    w = sqrt(t * u + 1)
    z = (1 + camera.alpha * (w - 1)) / v
    z_ = (camera.alpha * camera.beta) * ((t * v / w) + (2 * camera.alpha) * (camera.alpha * (w - 1) + 1)) / v
    return y * z, (W .* (camera.f' / z)) * (I - (y * y') * z_ / (z + n * z_))
end

function makeeucamera(coeffs, imhalfsz)
    # Convert the lens distortion 

end