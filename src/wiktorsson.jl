# Wiktorsson method
# First published by Wiktorsson, 2001

struct Wiktorsson <: AbstractIteratedIntegralAlgorithm end

convorder(::Wiktorsson) = 1//1
errcoeff(m, h, ::Wiktorsson, ::MaxL2) = √(5*m)*h/(√12*π)
norv(m, n, ::Wiktorsson) = 2*m*n+(m^2-m)÷2


"""
    levyarea(W, n, alg::Wiktorsson)

Simulates an approximation of the iterated Itô-integrals ``\\int_0^1W_s\\otimes dW_s``
of the given ``m``-dimensional increment of a Wiener process with step size 1.
The parameter ``n`` specifies the number of terms in the approximation and thus determines the accuracy.
This is an efficient implementation of the algorithm proposed in [Wiktorsson, 2001](@ref wiktorsson2001).
It is based on the Fourier method from Milstein but incorporates an additional tail sum approximation.
The algorithm needs approximately ``2\\cdot m^2+2\\cdot m\\cdot n+m`` Float64's.
The time complexity is ``\\mathcal{O}(m^2\\cdot n)``.
"""
function levyarea(W::AbstractVector{T}, n::Integer, alg::Wiktorsson) where {T<:AbstractFloat}
    # Preallocate
    m = length(W)
    A = similar(W,m,m) # allocates m*m Floats
    G = similar(W,m,m) # allocates m*m Floats
    # 1. Simulate Xₖ and Yₖ and approximate stochastic area integral
    Y = randn(m,n) # allocates m*n Floats
    Y .= (Y .- √(2).*W) ./ (1:n)'
    mul!(A,Y,randn(n,m)) # allocates m*n Floats
    # 2.a Simulate Gₙ
    a = √(2*trigamma(n+1))
    for j=1:m
        @inbounds G[j,j] = 0.0
        for i=j+1:m
            g = a * randn()
            @inbounds G[i,j] = g
            @inbounds G[j,i] = -g
            @inbounds A[i,j] += g
        end
    end
    # 2.b and add the tail-sum approximation
    A .+= inv(1+√(1+W'*W)) .* (G*W) .* W' # allocates m Floats
    # 3. Calculate the iterated integrals
    G .= inv(2pi).*(A .- A') # reuse G to save allocations
    return G
end