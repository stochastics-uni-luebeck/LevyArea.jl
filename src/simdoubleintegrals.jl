"""
    simdoubleintegrals(W::AbstractVector{AbstractFloat64}, n::Integer)

Simulates an approximation of all one-time iterated Itô-integrals
of the given Brownian motions with step size 1.
The algorithm is an adaptation of [^Wiktorsson2001].
The algorithm needs approximately ``2\\cdot m^2+2\\cdot m\\cdot n+m`` Float64's.
The time complexity should be around ``\\mathcal{O}(m^2\\cdot n)``.

Input:  W   the increments of m Brownian motions, where `m = length(W)`
        n   number of terms in the approximation of the stochastic area integral
Output: I[i,j] is an approximation of ``\\int_0^1W_i(s)dW_j(s)``

[^Wiktorsson2001]: *"Joint characteristic function and simultaneous simulation
    of iterated Itô integrals for multiple independent Brownian motions."*
    The Annals of Applied Probability 11.2: 470-487.

"""
function simdoubleintegrals_n(W::AbstractVector{<:AbstractFloat}, n::Integer;
    rng = GLOBAL_RNG)
    # Preallocate
    m::Integer = length(W)
    A = similar(W,m,m) # allocates m*m Floats
    G = similar(W,m,m) # allocates m*m Floats
    # 1. Simulate Brownian motion and determine n
    #    These are given as input arguments.
    # 2. Simulate Xₖ and Yₖ and approximate stochastic area integral
    Y = randn(rng,m,n) # allocates m*n Floats
    Y .= (Y .+ √(2).*W) ./ (1:n)'
    mul!(A,Y,randn(rng,n,m)) # allocates m*n Floats
    # 3.a Simulate Gₙ
    a = √(2*trigamma(n+1))
    for j=1:m
        @simd for i=1:j-1
            @inbounds G[i,j] = - @inbounds G[j,i]
        end
        @inbounds G[j,j] = 0.0
        @simd for i=j+1:m
            @inbounds G[i,j] = a * randn()
            @inbounds A[i,j] += @inbounds G[i,j]
        end
    end
    # 3.b and add the tail-sum approximation
    A .+= inv(1+√(1+W'*W)) .* W .* (W'*G) # allocates m Floats
    G .= (A .- A').*inv(2π) # reuse G to save allocations
    # 4. Calculate the iterated integrals
    BLAS.ger!(0.5,W,W,G) # 0.5*W*W' + G
    @inbounds for i=1:m # G-0.5I
        G[i,i] -= 0.5
    end
    G
end

"""
    simdoubleintegrals(W::AbstractVector{AbstractFloat64}, h::Real, C::Real=1.0)
    simdoubleintegrals(W::Real, h::Real=1.0)

Simulates an approximation of the iterated Itô-integrals
``\\int_0^h\\int_0^sdW_i(t)dW_j(s)`` for all pairs ``1\\le i, j \\le m``
of the given m-dimensional Brownian motion with step size h.
For the used algorithm see [^Wiktorsson2001].

In the case of a scalar Brownian motion the integral can be explicitly
calculated as ``\\int_0^h\\int_0^sdW(t)dW(s) = \\frac{1}{2}W(h)^2 - \\frac{1}{2}h``.

"""
function simdoubleintegrals(W::AbstractVector{<:AbstractFloat}, h::Real, C::Real=1.0; kwargs...)
    m::Integer = length(W)
    n::Int64 = ceil(Int64, √( m*(m-1)*(m+4*(W'*W)/h)/(C*h*24*π^2) ))
    Iᵢⱼ = h * simdoubleintegrals_n(W/√h, n; kwargs...)
    return Iᵢⱼ
end

simdoubleintegrals(W::Real, h::Real=1.0) = 0.5*W^2 - 0.5*h

####################################################

# struct DoubleIntegralCache
#     m::Int # dimension of Brownian motion
#     h::Real # stepsize
#     n::Int # number of terms in the approximation
#     a::Real
#     X::Matrix{Float64}
#     Y::Matrix{Float64}
#     A::Matrix{Float64}
#     G::Matrix{Float64}
#     function DoubleIntegralCache(m::Int,h::Real,C::Real)
#         n = ceil(Int, √( 5*m^2*(m-1)/(C*h*24*π^2) ))
#         a = √(2*trigamma(n+1))
#         X = Matrix{Float64}(undef,n,m)
#         Y = Matrix{Float64}(undef,m,n)
#         A = Matrix{Float64}(undef,m,m)
#         G = Matrix{Float64}(undef,m,m)
#         new(m,h,n,a,X,Y,A,G)
#     end
# end
# function simdoubleintegrals(W::AbstractVector{<:AbstractFloat},
#                             cache::DoubleIntegralCache;
#                             rng = GLOBAL_RNG)
#     # 1. Simulate Brownian motion and determine n
#     #    These are given as input arguments.
#     # 2. Simulate Xₖ and Yₖ and approximate stochastic area integral
#     randn!(rng,cache.X)
#     randn!(rng,cache.Y)
#     cache.Y .= (cache.Y .+ √(2) .* W)
#     cache.Y .= cache.Y ./ (1:cache.n)'
#     mul!(cache.A, cache.Y, cache.X)
#     # 3.a Simulate Gₙ
#     for j=1:cache.m
#         for i=1:j
#             @inbounds cache.G[i,j] = 0.0
#         end
#         for i=j+1:cache.m
#             @inbounds cache.G[i,j] = cache.a * randn()
#         end
#     end
#     # 3.b and add the tail-sum approximation
#     cache.A .+= 1 ./ (1+√(1+W'*W)) .* W.* (W'*cache.G .- W'*cache.G') .+ cache.G
#     cache.G .= cache.A .- cache.A' # reuse G to save allocations
#     # 4. Calculate the iterated integrals
#     cache.A .= (W*W'-I)./2 .+ cache.G./2π
#     cache.A
# end
#
# W = randn(5); C = DoubleIntegralCache(5,1,1);
# W100 = randn(100); C100 = DoubleIntegralCache(100,1,1);
# W1000 = randn(1000); C1000 = DoubleIntegralCache(1000,1,1);
#
# @benchmark simdoubleintegrals($W1000,4593)
# @benchmark simdoubleintegrals($W1000,$C1000)
