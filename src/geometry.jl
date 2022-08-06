export rodrigues
using StaticArrays

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
