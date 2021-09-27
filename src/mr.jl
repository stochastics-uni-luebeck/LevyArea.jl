# Mrongowius-Rößler method
# First published by Mrongowius & Rößler, 2021

struct MR{T<:AbstractFloat} <: AbstractIteratedIntegralAlgorithm
    m::Int
    n::Int
    X::Matrix{T}
    Y::Matrix{T}
    A::Matrix{T}
    MR{T}(m, n) where {T<:AbstractFloat} = new{T}(m,n,Matrix{T}(undef,m,n),Matrix{T}(undef,n,m),Matrix{T}(undef,m,m))
end

convorder(::MR) = 1//1
errcoeff(m, h, ::MR, ::MaxL2) = √m*h/(√12*π)
norv(m, n, ::MR) = 2*m*n+(m^2+m)÷2


"""
    levyarea(W, n, alg::MR)

Simulates an approximation of the iterated Itô-integrals ``\\int_0^1W_s\\otimes dW_s``
of the given ``m``-dimensional increment of a Wiener process with step size 1.
The parameter ``n`` specifies the number of terms in the approximation and thus determines the accuracy.
This is an efficient implementation of the algorithm proposed in [Mrongowius & Rößler, 2021](@ref mr2021).
It is based on the Fourier method from Milstein but incorporates an improved tail sum approximation.
The algorithm needs approximately ``m^2+2\\cdot m\\cdot n`` Float's 
and ``1/2m^2+2\\cdot m\\cdot n + 1/2m`` random numbers.
The time complexity is ``\\mathcal{O}(m^2\\cdot n)``.
"""
function levyarea(W::AbstractVector{T}, alg::MR{T}) where {T<:AbstractFloat}
    rng = default_rng()
    m = length(W)
    @assert m == alg.m
    n = alg.n
    X,Y,A = alg.X, alg.Y, alg.A
    # 1. Simulate Xₖ and Yₖ and approximate stochastic area integral
    randn!(rng, X)
    randn!(rng, Y)
    Y .= (Y .- √(T(2)).*W') ./ (1:n)
    mul!(A,X,Y)
    
    # 2. Add first, simple rest approximation (a₀)
    Ψ = randn!(rng, view(X, :, 1))
    a = T(sqrt(2*trigamma(n+1)))
    A .+= a .* W .* Ψ'

    # 3. Add improved tail-sum approximation and antisymmetrize (A.-=A')
    for i = 1:m
        @inbounds A[i,i] = zero(T)
        for j = (i+1):m
            @inbounds A[i,j] = (A[i,j] + a*randn(rng, T) - A[j,i])/(2*T(π))
            @inbounds A[j,i] = -A[i,j]
        end
    end

    return A
end
