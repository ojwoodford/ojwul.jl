export NoDistortionCamera, ExtendedUnifiedCamera
export ideal2image, image2ideal
abstract type AbstractCamera end
using StaticArrays

struct NoDistortionCamera{T<:Number} <: AbstractCamera
    fx::T
    fy::T
    cx::T
    cy::T
end

function ideal2image(camera::NoDistortionCamera, x)
    return SVector(x[1] * camera.fx + camera.cx, x[2] * camera.fy + camera.cy)
end

function image2ideal(camera::NoDistortionCamera, x)
    return SVector((x[1] - camera.cx) / camera.fx, (x[2] - camera.cy) / camera.fy)
end


struct ExtendedUnifiedCamera{T<:Number} <: AbstractCamera
    fx::T
    fy::T
    cx::T
    cy::T
    alpha::T
    beta::T
end

function ideal2image(camera::ExtendedUnifiedCamera, x)
    z = 1 / (1 + camera.alpha * (sqrt(camera.beta * (x[1] ^ 2 + x[2] ^ 2) + 1) - 1))
    return SVector(x[1] * z * camera.fx + camera.cx, x[2] * z * camera.fy + camera.cy)
end

function image2ideal(camera::ExtendedUnifiedCamera, x)
    y = SVector((x[1] - camera.cx) / camera.fx, (x[2] - camera.cy) / camera.fy)
    z = camera.beta * (y[1] ^ 2 + y[2] ^ 2)
    z = (1 + camera.alpha * (sqrt(1 + z * (1 - 2 * camera.alpha)) - 1)) / (1 - z * (camera.alpha ^ 2))
    return y * z
end
