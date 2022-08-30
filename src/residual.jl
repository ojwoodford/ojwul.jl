export cost, costresjac!

function cost(residual::AbstractResidual, vars::Vector{<:AbstractVariable})
    # Dispatch to the correct residual function with the correct block of variables
    r = computeresidual(residual, vars[varindices(residual)]...)
    # Compute the  robustified cost
    return robustify(r' * r, robustkernel(residual))
end

function costresjac!(res, jac, residual::AbstractResidual, vars::Vector{<:AbstractVariable})
    # Get the variables, compute the offsets
    
    # Convert to autodiff variables
    
    # Dispatch to the cost function
    r = cost(residual, vars...)

    # Extract the residuals and Jacobian from the autodiff output
    jac .= grad(r)
    res .= value(r)

    # Compute the cost and scale the residual and Jacobian 
    c, w = robustify(res' * res, robustkernel(residual))
    if w != 1
        res .*= w
        jac .*= w
    end
    return c
end