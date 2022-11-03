import ForwardDiff
export cost, costgradhess!

function cost(residuals::Vector{<:AbstractResidual}, vars::Vector{<:AbstractVariable})
    # Compute the total cost of all residuals
    c = 0.
    foreach(res -> c += cost(res, vars), residuals)
    return c
end

function cost(residual::AbstractResidual, vars::Vector{<:AbstractVariable})
    # Dispatch to the correct residual function with the correct block of variables
    r = computeresidual(residual, vars[varindices(residual)]...)
    # Compute the  robustified cost
    return robustify(robustkernel(residual), r' * r)[1]
end

function computejacobian(x, residual::AbstractResidual, vars, blockindex)
    # Update the variables we wish to differentiate
    start = 0
    vars = SizedArray(vars)
    for ind = range(1, length(vars))
        if blockindex[ind] != 0
            vars[ind] = update(vars[ind], x, start)
            start += nvars(vars[ind])
        end
    end
    # Compute the residual
    return computeresidual(residual, vars...)
end

function gradhesshelper!(grad, hess, fun, res, w, blockoffsets, zeros)
    # Compute the Jacobian
    jac = ForwardDiff.jacobian(fun, zeros)

    # IRLS weighting of the residual and Jacobian 
    if w != 1
        res = res * w
        jac = jac * w
    end

    # Update the blocks in the problem
    grad[blockoffsets] += jac' * res
    hess[blockoffsets, blockoffsets] += jac' * jac
    return nothing
end

function gradhesshelperfixed!(grad, hess, fun, res, w, blockoffsets, ::Val{N}) where N
    gradhesshelper!(grad, hess, fun, res, w, blockoffsets, zeros(SVector{N, eltype(res)}))
end

function costgradhess!(grad, hess, residual::AbstractResidual, vars::Vector{<:AbstractVariable}, blockindex::Vector{Int})
    # Get the variables
    varind = varindices(residual)
    vars = vars[varind]

    # Compute the robust residual
    res = computeresidual(residual, vars...)
    c, w = robustify(robustkernel(residual), res' * res)

    # Compute the number of variables and the block offsets
    blockindex = blockindex[varind]
    numvars = 0
    blockoffsets = Vector{Int}()
    for ind = range(1, length(varind))
        if blockindex[ind] != 0
            N = nvars(vars[ind])
            numvars += nvars(vars[ind])
        end
    end
    if numvars == 0
        # Bail out early if no variables
        return c
    end
    
    # Compute the hessian and gradient
    if numvars <= 32
        # Use a fixed array size
        valuedispatch1to32(v -> gradhesshelperfixed!(grad, hess, x -> computejacobian(x, residual, vars, blockindex), res, w, blockoffsets, v), numvars)
    else
        gradhesshelper!(grad, hess, x -> computejacobian(x, residual, vars, blockindex), res, w, blockoffsets, zeros(eltype(res), numvars))
    end

    # Return the cost
    return c
end
