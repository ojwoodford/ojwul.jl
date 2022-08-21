export NoDistortionCamera, ExtendedUnifiedCamera
export ideal2image, image2ideal, pixel2image, image2pixel
abstract type AbstractCamera end
using StaticArrays


function image2pixel(imsz, x)
    return SVector((x[1] + 0.5) * imsz[2], x[2] * imsz[2] + imsz[1] * 0.5)
end

function pixel2image(imsz, x)
    return SVector(x[1] / imsz[2] - 0.5, (x[2] - imsz[1] * 0.5) / imsz[2])
end

function pixel2image(imsz, x, W)
    return (SVector(x[1] / imsz[2] - 0.5, (x[2] - imsz[1] * 0.5) / imsz[2]), W ./ imsz[2])
end


struct SimpleCamera{T<:Number} <: AbstractCamera
    f::T
end

function ideal2image(camera::SimpleCamera, x)
    return x .* camera.f
end

function image2ideal(camera::SimpleCamera, x)
    return x ./ camera.f
end

function image2ideal(camera::SimpleCamera, x, W)
    return (x ./ camera.f, W ./ camera.f)
end


struct NoDistortionCamera{T<:Number} <: AbstractCamera
    f::SVector{2, T}
    c::SVector{2, T}
end

function ideal2image(camera::NoDistortionCamera, x)
    return x .* camera.f .+ camera.c
end

function image2ideal(camera::NoDistortionCamera, x)
    return (x .- camera.c) ./ camera.f
end


struct ExtendedUnifiedCamera{T<:Number} <: AbstractCamera
    f::SVector{2, T}
    c::SVector{2, T}
    alpha::T
    beta::T
end

function ideal2image(camera::ExtendedUnifiedCamera, x)
    z = 1 / (1 + camera.alpha * (sqrt(camera.beta * (x[1] ^ 2 + x[2] ^ 2) + 1) - 1))
    return (x .* z) .* camera.f .+ camera.c
end

function image2ideal(camera::ExtendedUnifiedCamera, x)
    y = (x .- camera.c) ./ camera.f
    z = camera.beta * (y[1] ^ 2 + y[2] ^ 2)
    z = (1 + camera.alpha * (sqrt(1 + z * (1 - 2 * camera.alpha)) - 1)) / (1 - z * (camera.alpha ^ 2))
    return y .* z
end
