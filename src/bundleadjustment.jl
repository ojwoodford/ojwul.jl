using VisualGeometryDatasets, ojwul
export BALImage, BALResidual
export getvars, computeresidual, robustkernel, nvars, transform, makeBALproblem

struct BALImage{T<:Real} <: AbstractVariable
    pose::EffPose3D{T}
    camera::BALCamera{T}
end
nvars(::BALImage) = 9
function update(var::BALImage, updatevec, start=0)
    return BALImage(update(var.pose, updatevec, start),
                    update(var.camera, updatevec, start+6))
end
function transform(im::BALImage, X::Point3D)
    return ideal2image(im.camera, -project(im.pose * X))
end
function makeBALImage(rx::T, ry::T, rz::T, tx::T, ty::T, tz::T, f::T, k1::T, k2::T) where T<:Real
    R = Rotation3DL(rx, ry, rz)
    return BALImage{T}(EffPose3D(R, Point3D(R.m' * -SVector(tx, ty, tz))), BALCamera(f, k1, k2))
end
makeBALImage(v) = makeBALImage(v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], v[9])

struct BALResidual{T<:Real} <: AbstractResidual
    measurement::SVector{2, T}
    varind::SVector{2, Int}
end
BALResidual(m, v) = BALResidual(SVector{2}(m[1], m[2]), SVector{2, Int}(v[1], v[2]))
nvars(::BALResidual) = 2 # Residual depends on 2 variables
reslen(::BALResidual) = 2 # The residual is a vector of length 2
function eltype(::BALResidual{T}) where T
    return T
end
function getvars(res::BALResidual{T}, vars::Vector{<:AbstractVariable}) where T
    return vars[res.varind[1]]::BALImage{T}, vars[res.varind[2]]::Point3D{T}
end
function computeresidual(res::BALResidual, im::BALImage, X::Point3D)
    return transform(im, X) - res.measurement
end
const balrobustifier = HuberKernel(2., 4., 1., 1.)
function robustkernel(::BALResidual)
    return balrobustifier
end


function makeBALproblem(name)
    # Load the data
    data = loadbaldataset(name)

    # Construct the cameras and landmarks
    vars = Vector{AbstractVariable}()
    # Add the cameras
    for col = eachcol(data.cameras)
        push!(vars, makeBALImage(col))
    end
    numcameras = length(vars)
    # Add the landmarks
    for col = eachcol(data.landmarks)
        push!(vars, Point3D(col[1], col[2], col[3]))
    end

    # Construct the residuals
    residuals = Vector{BALResidual{Float64}}()
    for ind = 1:size(data.measurementindices, 2)
        push!(residuals, BALResidual(data.measurements[:,ind], SVector{2, Int}(data.measurementindices[1,ind], data.measurementindices[2,ind] + numcameras)))
    end

    # Return the optimization problem
    return (vars, residuals)
end
