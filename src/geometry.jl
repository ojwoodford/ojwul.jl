export rodrigues
using StaticArrays
using LinearAlgebra

function rodrigues(x::T, y::T, z::T) where T<:AbstractFloat
theta2 = x * x + y * y + z * z
cosf = T(0.5);
sinc = T(1);
if theta2 > T(2.23e-16)
    theta = sqrt(theta2)
    cosf = (T(1) - cos(theta)) / theta2;
    sinc = sin(theta) / theta;
end
out = MMatrix{3, 3, T}(undef)
# Diagonals
out[1,1] = (x * x - theta2) * cosf + 1;
out[2,2] = (y * y - theta2) * cosf + 1;
out[3,3] = (z * z - theta2) * cosf + 1;
# Off diagonals
a = x * y * cosf
b = sinc * z
out[2,1] = a + b;
out[1,2] = a - b;
a = x * z * cosf
b = sinc * y
out[3,1] = a - b;
out[1,3] = a + b;
a = y * z * cosf
b = sinc * x
out[2,3] = a - b;
out[3,2] = a + b;
return out
end