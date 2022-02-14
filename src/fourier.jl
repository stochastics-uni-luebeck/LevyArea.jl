# Basic Fourier method without any kind of rest approximation

struct Fourier <: AbstractIteratedIntegralAlgorithm end

convorder(::Fourier) = 1//2
errcoeff(m, h, ::Fourier, ::MaxL2) = h*√3/(√2*π)
norv(m, n, ::Fourier) = 2*m*n


"""
    levyarea(W, n, alg::Fourier)

Simulates an approximation of the iterated Itô-integrals ``\\int_0^1W_s\\otimes dW_s``
of the given ``m``-dimensional increment of a Wiener process with step size 1.
The parameter ``n`` specifies the number of terms in the approximation and thus determines the accuracy.
This algorithm is based on a Fourier expansion of the Wiener process.
The algorithm needs approximately ``m^2+2\\cdot m\\cdot n`` Float's
and ``2\\cdot m\\cdot n`` random numbers.
The time complexity is ``\\mathcal{O}(m^2\\cdot n)``.
"""
function levyarea(W::AbstractVector{T}, n::Integer, alg::Fourier) where {T}
    rng = default_rng()
    m = length(W)
    X = randn(rng, T, n, m) # allocates m*n Floats
    Y = randn(rng, T, m, n) # allocates m*n Floats
    Y .= (Y .- √(T(2)).*W) ./ (1:n)'
    A = Y * X # allocates m*m Floats
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
