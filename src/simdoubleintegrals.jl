"""
    simdoubleintegrals(W::Vector{Float64}, n::Integer)

Simulates an approximation of all one-time iterated Itô-integrals
of the given Brownian motions with step size 1.
The algorithm is taken from [^Wiktorsson2001].

Input:  W   the increments of m Brownian motions, where `m = length(W)`
        n   number of terms in the approximation of the stochastic area integral
Output: I[i,j] is an approximation of ``\\int_0^hW_i(s)dW_j(s)``

### References
[^Wiktorsson2001]: *"Joint characteristic function and simultaneous simulation
    of iterated Itô integrals for multiple independent Brownian motions."*
    The Annals of Applied Probability 11.2: 470-487.
"""
function simdoubleintegrals(W::Vector{Float64}, n::Integer)
    # Preparation
    m::Integer = length(W)
    M::Integer = m*(m-1) ÷ 2
    Pₘ = sparse(1:m^2, cld.(1:m^2, m).+m.*mod.(0:(m^2-1), m), 1)
    Kₘ = sparse(1:M, [i+(j-1)*m for i=1:m, j=1:m if j<i], 1, M, m^2)
    KIP = Kₘ*(I-Pₘ)

    # 1. Simulate Brownian motion and determine n
    # These are given as input arguments.
    # 2.a Simulate Xₖ and Yₖ
    X = randn(m,n)
    YW = randn(m,n) .+ √(2).*W
    # 2.b Approximate stochastic area integral
    Aₙ = X*Diagonal(1 ./ (1:n))*YW'
    Aₙ = Aₙ'-Aₙ
    # 3. Simulate Gₙ and add the tail-sum approximation
    temp = KIP * kron(sparse(1.0I,m,m),W)
    Σ = 2.0 * temp * temp' + 2.0I
    sqrtΣ = (Σ + (2*√(1+W'*W))*I)
    lmul!(1 / (√2 + √(2+2*W'*W)), sqrtΣ)
    Aₙ .= Aₙ .+ √trigamma(n+1) .* reshape( KIP'*sqrtΣ*randn(M) ,m,m)
    # 4. Calculate the iterated integrals
    Iᵢⱼ = (W*W'-I) ./2 .+ 1/2π .* Aₙ
    return Iᵢⱼ
end

function simdoubleintegrals(W::Vector{Float64}, h::Real, C::Real)
    lmul!(1/√h, W)
    n::Int64 = ceil(Int64, √( m*(m-1)*(m+4*(W'*W))/(C*h*24*π^2) ))
    Iᵢⱼ = simdoubleintegrals(W, n)
    lmul!(√h, W)
    lmul!(h, Iᵢⱼ)
    return Iᵢⱼ
end
