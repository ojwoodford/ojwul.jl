import ForwardDiff, Zygote, ojwul.valuedispatch
export cost, costgradhess!

function cost(residuals, vars::Vector{<:AbstractVariable})
    # Compute the total cost of all residuals
    c = 0.
    for res in residuals
        c += cost(res, vars)
    end
    return c
end

function cost(residual::Residual, vars::Vector{<:AbstractVariable}) where Residual <: AbstractResidual
    # Get the variables required to compute the residual
    v = getvars(residual, vars)
    # Compute the residual
    r = computeresidual(residual, v...)
    # Compute the robustified cost
    return robustify(robustkernel(residual), r' * r)[1]
end

# Count the numbers of variables in the inputs
function countvars(vars::Tuple, ::Val{varflags}) where varflags
    flags = varflags
    count = 0
    for v in vars
        count += nvars(v) * (flags & 1)
        flags >>= 1
    end
    return count
end

# Generate the updated variables
function updatevars(vars::Tuple, ::Val{varflags}, advar) where varflags
    return ntuple(i -> ((varflags >> (i - 1)) & 1) != 0 ? update(vars[i], advar, countvars(vars[1:i-1], Val(varflags))) : vars[i], length(vars))
end

# Compute the offsets of the variables
function computeoffsets(vars::Tuple, ::Val{varflags}, blockind) where varflags
    return vcat(ntuple(i -> SR(1, nvars(vars[i]) * ((varflags >> (i - 1)) & 1)) .+ (blockind[i] - 1), length(vars))...)
end

function gradhesshelperfd!(grad, hess, residual::Residual, vars::Vector{<:AbstractVariable}, blockind, ::Val{varflags}) where {varflags, Residual <: AbstractResidual}
    # Get the variables
    v = getvars(residual, vars)

    # Compute the residual
    res = computeresidual(residual, v...)

    # Compute the robustified cost and the IRLS weight
    c, w = robustify(robustkernel(residual), res' * res)

    # Bail early if the weight is zero
    if w == 0
        return c
    end

    # Compute the Jacobian
    jac = ForwardDiff.jacobian(z -> computeresidual(residual, updatevars(v, Val(varflags), z)...), zeros(SVector{countvars(v, Val(varflags)), eltype(res)}))

    # IRLS weighting of the residual and Jacobian 
    if w != 1
        res = res * w
        jac = jac * w
    end

    # Update the blocks in the problem
    blockoffsets = computeoffsets(v, Val(varflags), blockind)
    grad[blockoffsets] += jac' * res
    hess[blockoffsets, blockoffsets] += jac' * jac

    # Return the cost
    return c
end

function gradhesshelperzy!(grad, hess, residual::Residual, vars::Vector{<:AbstractVariable}, blockind, ::Val{varflags}) where {varflags, Residual <: AbstractResidual}
    # Get the variables
    v = getvars(residual, vars)

    # Compute the Jacobian
    res, jac = Zygote.forward_jacobian(z -> computeresidual(residual, updatevars(v, Val(varflags), z)...), zeros(SVector{countvars(v, Val(varflags)), Float64}))

    # Compute the robustified cost and the IRLS weight
    c, w = robustify(robustkernel(residual), res' * res)

    # Bail early if the weight is zero
    if w == 0
        return c
    end

    # IRLS weighting of the residual and Jacobian 
    if w != 1
        res = res * w
        jac = jac * w
    end

    # Update the blocks in the problem
    blockoffsets = computeoffsets(v, Val(varflags), blockind)
    grad[blockoffsets] += jac * res
    hess[blockoffsets, blockoffsets] += jac * jac'

    # Return the cost
    return c
end

function costgradhess!(grad, hess, residuals, vars::Vector{<:AbstractVariable}, blockindex::Vector{Int})
    # Go over all resdiduals, updating the gradient & hessian, and aggregating the cost 
    c = 0.
    for res in residuals
        c += costgradhess!(grad, hess, res, vars, blockindex)
    end
    return c
end

function costgradhess!(grad, hess, residual::Residual, vars::Vector{<:AbstractVariable}, blockindex::Vector{Int}) where Residual <: AbstractResidual
    # Get the bitset for the input variables, as an integer
    blockind = blockindex[residual.varind]
    varflags = foldl((x, y) -> (x << 1) + (y != 0), blockind, init=0)
    # If there are no variables, just return the cost
    if varflags == 0
        return cost(residual, vars)
    end
    # Dispatch gradient computation based on the varflags, and return the cost
    #return gradhesshelperfd!(grad, hess, residual, vars, blockind, Val(varflags))
    return valuedispatch(Val(1), Val((2^nvars(residual))-1), v -> gradhesshelperfd!(grad, hess, residual, vars, blockind, v), varflags)
end
