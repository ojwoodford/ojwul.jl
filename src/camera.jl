export NoDistortionCamera, ExtendedUnifiedCamera
export ideal2image, image2ideal, pixel2image, image2pixel
abstract type AbstractCamera end
using StaticArrays


function image2pixel(halfimsz, x)
    return x .* halfimsz[1] .+ halfimsz
end

function pixel2image(halfimsz, x)
    return (x .- halfimsz) ./ halfimsz[1]
end

function pixel2image(halfimsz, x, W)
    return ((x .- halfimsz) ./ halfimsz[1], W ./ halfimsz[1])
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
    z = 1 ./ (1 .+ camera.alpha .* (sqrt.(camera.beta .* sum(abs2, x, dims=1) .+ 1) .- 1))
    return (x .* z) .* camera.f .+ camera.c
end

function image2ideal(camera::ExtendedUnifiedCamera, x)
    y = (x .- camera.c) ./ camera.f
    z = camera.beta .* sum(abs2, y, dims=1)
    z = (1 .+ camera.alpha .* (sqrt.(1 .+ z .* (1 .- 2 .* camera.alpha)) .- 1)) ./ (1 .- z .* (camera.alpha ^ 2))
    return y .* z
end
