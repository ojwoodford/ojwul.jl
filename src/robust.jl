export robustify!

struct HuberKernel{T<:Real}
    width::T
    width_squared::T
end

HuberKernel(x) = HuberKernel(x, x*x)

function robustify!(cost::T, kernel::HuberKernel{T}) where T<:AbstractFloat
grad = ones(T, size(cost))
outliers = cost .> kernel.width_squared
sqrt_outliers = sqrt.(cost[outliers])
cost[outliers] .= sqrt_outliers .* (kernel.width * 2) .- kernel.width_squared
grad[outliers] .= kernel.width ./ sqrt_outliers
return grad
end