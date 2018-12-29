"""
    createBrownianPath(I::AbstractRange,m::Integer)

Generate an `m`-dimensional Brownian path at the timepoints specified by `I`.
The path always starts at ``0_m``.

"""
createBrownianPath(I::AbstractRange,m::Integer) = √step(I) * [zeros(1,m);randn(length(I)-1,m)]

"""
    restrict(I::AbstractRange,dW::AbstractArray,k::Integer=1)

Restrict a path of a Brownian motion given by a range of time points and the
corresponding Brownian increments to ``\\frac{1}{2^k}`` as many time points.
The lengths of `I` and the first dimension of `dW` must match and be of the
form ``2^n+1`` for some ``n\\ge k``.

"""
function restrict(I::AbstractRange,dW::AbstractArray,k::Integer=1)
    l = length(I)
    l ≠ size(dW,1) && error("Lengths do not match!")
    !ispow2(l-1) && error("Length must be 2^n+1!")
    l < 2^k+1 && error("Length must be at least 2^k+1!")

    new_len = (l-1) ÷ (2^k) + 1
    I2 = I[1:2^k:end]
    dW2 = similar(dW,new_len,size(dW,2))
    dW2[1,:] .= dW[1,:]
    for i = 2:new_len
        dW2[i,:] .= 0.0
        for j = (i-2)*2^k+2:(i-1)*2^k+1
            dW2[i,:] .+= dW[j,:]
        end
    end

    (I2,dW2)
end
