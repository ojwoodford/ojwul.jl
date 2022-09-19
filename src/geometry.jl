export rodrigues, proj, homg, epipolarerror, proj2orthonormal
export Rotation3D, Point3D, Pose3D, UnitPose3D
using StaticArrays, LinearAlgebra
import ojwul.AbstractVariable, ojwul.SR

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

function proj2orthonormal(M)
    s = svd(M);
    return s.U * s.V';
end

function proj(x)
    return x[1:end-1,:] ./ x[end,:]
end

function proj(x::StaticVector)
    return x[SR(1, end-1)] ./ x[end]
end

function proj(x::StaticArray)
    return x[SR(1, end-1),:] ./ x[end,:]
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


struct Point3D{T<:Real} <: AbstractVariable
    v::SVector{3, T}
end
Point3D(x, y, z) = Point3D(SVector{3}(x, y, z))
Point3D() = Point3D(SVector{3}(0., 0., 0.))
function nvars(var::Point3D)
    return 3
end
function update(var::Point3D, updatevec)
    return Point3D(var.v + updatevec)
end
function update(var::Point3D, x, y, z)
    return Point3D(var.v + SVector(x, y, z))
end

function proj(x::Point3D)
    return SVector(x.v[1], x.v[2]) ./ x.v[3]
end


struct Rotation3D{T<:Real} <: AbstractVariable
    m::SMatrix{3, 3, T, 9}
end
Rotation3D(x, y, z) = Rotation3D(rodrigues(x, y, z))
Rotation3D() = Rotation3D(SMatrix{3, 3, Float64}(1., 0., 0., 0., 1., 0., 0., 0., 1.))
function nvars(var::Rotation3D)
    return 3
end
function update(var::Rotation3D, updatevec)
    return var * Rotation3D(updatevec[1], updatevec[2], updatevec[3])
end
function update(var::Rotation3D, x, y, z)
    return var * Rotation3D(x, y, z)
end
Base.:*(rota::Rotation3D, rotb::Rotation3D) = Rotation3D(rota.m * rotb.m)
Base.:*(rot::Rotation3D, point::Point3D) = Point3D(rot.m * point.v)


struct Pose3D{T<:Real} <: AbstractVariable
    rot::Rotation3D{T}
    trans::Point3D{T}
end
Pose3D(rx, ry, rz, tx, ty, tz) = Pose3D(Rotation3D(rx, ry, rz), Point3D(tx, ty, tz))
Pose3D() = Pose3D(Rotation3D(), Point3D())
function nvars(var::Pose3D)
    return 6
end
function update(var::Pose3D, updatevec)
    return Pose3D(update(var.rot, updatevec[1], updatevec[2] , updatevec[3]), 
                  update(var.trans, updatevec[SR(4, 6)]))
end
function inverse(var::Pose3D)
    return Pose3D(var.rot', var.rot' * -var.trans)
end
Base.:*(pose::Pose3D, point::Point3D) = Point3D(pose.rot.m * point.v + pose.trans.v)


struct UnitPose3D{T<:Real} <: AbstractVariable
    rot::Rotation3D{T}
    trans::Rotation3D{T}
end
UnitPose3D() = Pose3D(Rotation3D(), Rotation3D())
UnitPose3D((rx, ry, rz, tx, ty, tz)) = Pose3D(Rotation3D(rx, ry, rz), Rotation3D()) # Normalize translation and initialize y & z axes
function nvars(var::UnitPose3D)
    return 5
end
function update(var::UnitPose3D, updatevec)
    return UnitPose3D(update(var.rot, updatevec[1], updatevec[2], updatevec[3]), update(var.trans, 0, updatevec[4], updatevec[5]))
end
function inverse(var::UnitPose3D)
    return Pose3D(var.rot', var.rot' * -var.trans.m[:,1])
end
Base.:*(pose::UnitPose3D, point::Point3D) = Point3D(pose.rot.m * point.v + pose.trans.m[:,1])
