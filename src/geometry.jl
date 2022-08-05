export rodrigues
using StaticArrays
using LinearAlgebra

function rodrigues(x::T, y::T, z::T) where T<:AbstractFloat
theta2 = x * x + y * y + z * z
cosf = T(0.5)
sinc = T(1)
if theta2 > T(2.23e-16)
    theta = sqrt(theta2)
    cosf = (T(1) - cos(theta)) / theta2
    sinc = sin(theta) / theta
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