export rodrigues, proj, homg, epipolarerror
export Rotation3D, Point3D, Pose3D
using StaticArrays
using LinearAlgebra

function rodrigues(x::T, y::T, z::T) where T<:Number
    if x == 0 && y == 0 && z == 0
        # Short cut for derivatives at identity
        return SMatrix{3, 3, T}(T(1), z, -y, -z, T(1), x, y, -x, T(1))
    end
    theta2 = x * x + y * y + z * z
    cosf = T(0.5)
    sinc = T(1)
    if theta2 > T(2.23e-16)
        theta = sqrt(theta2)
        sinc, cosf = sincos(theta)
        cosf -= 1
        sinc /= theta
        cosf /= -theta2
    end
    a = x * y * cosf
    b = sinc * z
    c = x * z * cosf
    d = sinc * y
    e = y * z * cosf
    f = sinc * x
    return SMatrix{3, 3, T}((x * x - theta2) * cosf + 1, a + b, c - d,
                            a - b, (y * y - theta2) * cosf + 1, e + f,
                            c + d, e - f, (z * z - theta2) * cosf + 1)
end

function proj(x)
    return x[1:end-1,:] ./ x[end,:]
end

function proj(x::StaticVector)
    return x[SVector{end-1, Int}(1:end-1)] ./ x[end]
end

function proj(x::StaticArray)
    return x[SVector{end-1, Int}(1:end-1),:] ./ x[end,:]
end

function homg(x)
    return vcat(x, ones(typeof(x[1]), 1, size(x, 2)))
end

function homg(x::StaticVector)
    return vcat(x, 1)
end

function homg(x::StaticArray)
    return vcat(x, ones(SMatrix{1, Size(x)[2]}))
end

function epipolarerror(RX, T, x, W=Nothing)
    Tp = proj(T)
    xe = x - Tp
    RXp = proj(RX) - Tp
    RXp = SVector(RXp[2], -RXp[1])
    if W != Nothing
        RXp = W' \ RXp
        xe = W * xe
    end
    return dot(xe, normalize(RXp))
end


struct Point3D{T} <: FieldVector{3, T} <: AbstractVariable where T<:Real
    x::T
    y::T
    z::T
end
function ndims(var::Point3D)
    return 3
end
function update(var::Point3D, updatevec)
    return Point3D(var + updatevec)
end


struct Rotation3D{T<:Real} <: AbstractVariable
    m::SMatrix{3, 3, T}
end
function ndims(var::Rotation3D)
    return 3
end
function update(var::Rotation3D, updatevec)
    return Rotation3D(var.m * rodrigues(updatevec[1], updatevec[2], updatevec[3]))
end
function update(var::Rotation3D, x, y, z)
    return Rotation3D(var.m * rodrigues(x, y, z))
end


struct Pose3D{T<:Real} <: AbstractVariable
    rot::Rotation3D{T}
    trans::Point3D{T}
end
function ndims(var::Pose3D)
    return 6
end
function update(var::Pose3D, updatevec)
    return Pose3D(update(var.rot, updatevec[1:3]), update(var.trans, updatevec[4:6]))
end
function inverse(var::Pose3D)
    return Pose3D(var.rot', var.rot' * -var.trans)
end
Base.:*(pose::Pose3D, point::Point3D) = Point3D(pose.rot * point + pose.trans)


struct UnitPose3D{T<:Real} <: AbstractVariable
    rot::Rotation3D{T}
    trans::Rotation3D{T}
end
function ndims(var::UnitPose3D)
    return 5
end
function update(var::UnitPose3D, updatevec)
    return UnitPose3D(update(var.rot, updatevec[1:3]), update(var.trans, 0, updatevec[4], updatevec[5]))
end
function inverse(var::UnitPose3D)
    return Pose3D(var.rot', var.rot' * -var.trans.m[:,1])
end
Base.:*(pose::UnitPose3D, point::Point3D) = Point3D(pose.rot * point + pose.trans.m[:,1])
