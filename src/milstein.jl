# Fourier method
# First published by Milstein, 1994

struct Milstein <: AbstractIteratedIntegralAlgorithm end

convorder(::Milstein) = 1//2
errcoeff(m, h, ::Milstein, ::MaxL2) = h/(√2*π)
norv(m, n, ::Milstein) = 2*m*n+m


"""
    levyarea(W, n, alg::Milstein)

Simulates an approximation of the iterated Itô-integrals ``\\int_0^1W_s\\otimes dW_s``
of the given ``m``-dimensional increment of a Wiener process with step size 1.
The parameter ``n`` specifies the number of terms in the approximation and thus determines the accuracy.
This is an efficient implementation of the algorithm proposed in [Milstein, 1994](@ref milstein1994).
It is based on a Fourier expansion of the Wiener process.
The algorithm needs approximately ``m^2+2\\cdot m\\cdot n`` Float's.
The time complexity is ``\\mathcal{O}(m^2\\cdot n)``.
"""
function levyarea(W::AbstractVector{T}, n::Integer, alg::Milstein) where {T<:AbstractFloat}
    rng = default_rng()
    m = length(W)
    X = randn(rng, T, n, m) # allocates m*n Floats
    Y = randn(rng, T, m, n) # allocates m*n Floats
    Y .= (Y .- √(T(2)).*W) ./ (1:n)'
    A = Y * X # allocates m*m Floats
    # Add simple rest approximation (a₀)
    # M = randn(rng, T, m) # allocates m Floats
    M = randn!(rng, view(Y, :, 1))
    a = T(√(2*trigamma(n+1)))
    A .+= a .* W .* M'
    # Antisymmetrize
    # G = inv(2*T(π)).*(A .- A') # allocates m*m Floats
    for i = 1:m
        @inbounds A[i,i] = zero(T)
        for j = (i+1):m
            @inbounds A[i,j] = (A[i,j] - A[j,i])/(2*T(π))
            @inbounds A[j,i] = -A[i,j]
        end
    end
    return A
end
