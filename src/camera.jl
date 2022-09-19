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

struct EULensDistortion{T<:Number} <: AbstractVariable
    alpha::T
    beta::T
end
function nvars(var::EULensDistortion)
    return 2
end
function update(var::EULensDistortion, updatevec)
    alpha = var.alpha * exp(updatevec[1])
    alpha = alpha / (1 + (alpha - var.alpha))
    beta = (var.beta != 0 ? var.beta : eps(var.beta)) * exp(updatevec[2])
    return EULensDistortion(alpha, beta)
end

function ideal2distorted(lens::EULensDistortion, x)
    z = 1 ./ (1 .+ lens.alpha .* (sqrt.(lens.beta .* sum(abs2, x, dims=1) .+ 1) .- 1))
    return x .* z
end

function ideal2distorted(lens::EULensDistortion, x::StaticVector)
    z = 1 / (1 + lens.alpha * (sqrt(lens.beta * (x' * x) + 1) - 1))
    return x * z
end

function distorted2ideal(lens::EULensDistortion, x)
    z = lens.beta .* sum(abs2, x, dims=1)
    z = (1 .+ lens.alpha .* (sqrt.(1 .+ z .* (1 .- 2 .* lens.alpha)) .- 1)) ./ (1 .- z .* (lens.alpha ^ 2))
    return x .* z
end

function distorted2ideal(lens::EULensDistortion, x, W)
    t = 1 - 2 * lens.alpha
    n = x' * x
    u = lens.beta * n
    v = 1 - u * (lens.alpha ^ 2)
    w = sqrt(t * u + 1)
    z = (1 + lens.alpha * (w - 1)) / v
    z_ = (lens.alpha * lens.beta) * ((t * v / w) + (2 * lens.alpha) * (lens.alpha * (w - 1) + 1)) / v
    return x * z, W * (I / z - (x * x') * (z_ / (z * (z + n * z_))))
end

struct ExtendedUnifiedCamera{T<:Number} <: AbstractCamera
    sensor::NoDistortionCamera{T}
    lens::EULensDistortion{T}
end
ExtendedUnifiedCamera(f, c, a, b) = ExtendedUnifiedCamera(NoDistortionCamera(f, c), EULensDistortion(a, b))
function nvars(var::ExtendedUnifiedCamera)
    return 6
end
function update(var::ExtendedUnifiedCamera, updatevec)
    return ExtendedUnifiedCamera(update(var.sensor, updatevec[SR(1, 4)]),
                                 update(var.lens, updatevec[SR(5, 6)]))
end

function ideal2image(camera::ExtendedUnifiedCamera, x)
    return ideal2image(camera.sensor, ideal2distorted(camera.lens, x))
end

function image2ideal(camera::ExtendedUnifiedCamera, x)
    return distorted2ideal(camera.lens, image2ideal(camera.sensor, x))
end

function image2ideal(camera::ExtendedUnifiedCamera, x, W)
    y, W_ = image2ideal(camera.sensor, x, W)
    return distorted2ideal(camera.lens, y, W_)
end

struct LensDistortResidual{T<:Number} <: AbstractResidual
    rlinear::T
    rdistort::T
end
function computeresidual(residual::LensDistortResidual, lens::EULensDistortion)
    return rdistort - ideal2distorted(lens, rlinear)
end
function varindices(residual::LensDistortResidual)
    return 1
end

function makeeucamera(coeffs)
    # Compute the sensor values
    #halfimsz = 

    # Create an optimization problem to convert the lens distortion
    #residuals = 

    

end