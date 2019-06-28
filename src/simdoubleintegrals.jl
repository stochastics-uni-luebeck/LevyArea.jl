"""
    ito_correction!(I, h=1)

Applies the Itô-correction for iterated integrals to `I`.
This amounts to subtracting ``\\frac{1}{2}h`` from every element on the diagonal.

# Example
```jldoctest
julia> M = ones(5,5);

julia> IteratedIntegrals.ito_correction!(M, 0.5)

julia> M
5×5 Array{Float64,2}:
 0.75  1.0   1.0   1.0   1.0
 1.0   0.75  1.0   1.0   1.0
 1.0   1.0   0.75  1.0   1.0
 1.0   1.0   1.0   0.75  1.0
 1.0   1.0   1.0   1.0   0.75
```

"""
@inline function ito_correction!(I, h=1)
    @inbounds for i=axes(I,1) # hope we don't get non-square matrices
        I[i,i] -= h/2
    end
end


"""
    simdoubleintegrals(W, n, alg::Fourier)

Simulates an approximation of the iterated Itô-integrals ``\\int_0^1W_s\\otimes dW_s``
of the given ``m``-dimensional increment of a Wiener process with step size 1.
The parameter ``n`` specifies the number of terms in the approximation and thus determines the accuracy.
This is an efficient implementation of the algorithm proposed in [Milstein, 1994](@ref milstein1994).
It is based on a Fourier expansion of the Wiener process.
The algorithm needs approximately ``2\\cdot m^2+2\\cdot m\\cdot n+m`` Float64's.
The time complexity is ``\\mathcal{O}(m^2\\cdot n)``.
"""
function simdoubleintegrals(W::AbstractVector{<:AbstractFloat}, n::Integer, alg::Fourier)
    m::Integer = length(W)
    Y = randn(m,n) # allocates m*n Floats
    Y .= (Y .- √(2).*W) ./ (1:n)'
    A = Y*randn(n,m) # allocates m*n Floats
    M = randn(m) # allocates m Floats
    A .+= √(2*trigamma(n+1)) .* W.*M'
    G = 0.5.*W.*W' .+ inv(2pi).*(A .- A') # allocates m*m Floats
    return  G
end

"""
    simdoubleintegrals(W, n, alg::Wiktorsson)

Simulates an approximation of the iterated Itô-integrals ``\\int_0^1W_s\\otimes dW_s``
of the given ``m``-dimensional increment of a Wiener process with step size 1.
The parameter ``n`` specifies the number of terms in the approximation and thus determines the accuracy.
This is an efficient implementation of the algorithm proposed in [Wiktorsson, 2001](@ref wiktorsson2001).
It is based on the Fourier method from Milstein but incorporates an additional tail sum approximation.
The algorithm needs approximately ``2\\cdot m^2+2\\cdot m\\cdot n+m`` Float64's.
The time complexity is ``\\mathcal{O}(m^2\\cdot n)``.
"""
function simdoubleintegrals(W::AbstractVector{<:AbstractFloat}, n::Integer, alg::Wiktorsson)
    # Preallocate
    m::Integer = length(W)
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
    G .= 0.5.*W.*W' .+ inv(2pi).*(A .- A') # reuse G to save allocations
    return G
end

"""
    simdoubleintegrals(W::AbstractVector, h, eps; ito_correction=true, alg=Wiktorsson())

Simulates an approximation of the iterated stochastic integrals
``\\int_0^h\\int_0^sdW_i(t)dW_j(s)`` for all pairs ``1\\le i, j \\le m``
of the given ``m``-dimensional Brownian motion with step size h.

# Examples
```jldoctest
julia> h = 1/2;

julia> W = [1.0, 0.5]
2-element Array{Float64,1}:
 1.0
 0.5

julia> diag(simdoubleintegrals(W, h, h^(3/2))) ≈ 0.5*W.^2 .- 0.5h
true
```

    simdoubleintegrals(W::AbstractVector, q_12::AbstractVector, h, eps; ito_correction=true, alg=Wiktorsson())

Simulates an approximation of the iterated stochastic integrals for finite-dimensional approximations of
a Q-Wiener process with covariance matrix ``Q = Q^\\frac{1}{2}*Q^\\frac{1}{2}``.
Here `q_12` is a vector of the eigenvalues of ``Q^\\frac{1}{2}``; the square root of the covariance matrix.
Equivalently these are the square roots of the eigenvalues of ``Q``.
"""
function simdoubleintegrals(W::AbstractVector{<:AbstractFloat}, h::Real, eps::Real;
                            ito_correction=true,
                            alg::AbstractIteratedIntegralAlgorithm=Fourier())
    n = terms_needed(W, h, eps, alg)
    I = h * simdoubleintegrals(W/√h, n, alg)
    if ito_correction
        ito_correction!(I,h)
    end
    return I
end

function simdoubleintegrals(W::AbstractVector{<:AbstractFloat}, q_12::AbstractVector, h::Real, eps::Real;
                            ito_correction=true,
                            alg::AbstractIteratedIntegralAlgorithm=Fourier())
    n = terms_needed(W, q_12, h, eps, alg)
    I = simdoubleintegrals(W./q_12./√h, n, alg)
    if ito_correction
        ito_correction!(I)
    end
    return h.*q12'.*I.*q12 # scale correctly
end

"""
    simdoubleintegrals(W::Real, h::Real=1.0, eps::Real=0.0; ito_correction=true, kwargs...)

In the case of a scalar Brownian motion the integral can be explicitly
calculated as ``\\int_0^h\\int_0^sdW(t)dW(s) = \\frac{1}{2}W(h)^2 - \\frac{1}{2}h``.

The parameter `eps` (as well as all additional keyword arguments) has no effect but is available to provide the same interface as the multidimensional version.
"""
simdoubleintegrals(W::Real, h::Real=1.0, eps::Real=0.0; ito_correction=true, kwargs...) = ito_correction ? 0.5W^2 - 0.5h : 0.5W^2

"""
    terms_needed(W, stepsize, eps, alg)

Returns the number of terms in the approximating sum that is needed to ensure an ``L^2``-error of at most `eps`.
This depends on the current stepsize and the chosen algorithm.
For some algorithms this additionally depends on the realisation of the Wiener increment.

# Examples
```jldoctest
julia> h = 1/128;

julia> W = sqrt(h)*randn(10);

julia> IteratedIntegrals.terms_needed(W, h, h^(3/2), IteratedIntegrals.Fourier())
7
```

    terms_needed(W, q_12, stepsize, eps, alg)

Used for finite-dimensional approximations of a Q-Wiener process with covariance matrix
``Q = Q^\\frac{1}{2}*Q^\\frac{1}{2}``. Here `q_12` is a vector of the eigenvalues of ``Q^\\frac{1}{2}``;
the square root of the covariance matrix. Equivalently these are the square roots of the eigenvalues of ``Q``.
"""
function terms_needed end
terms_needed(W, stepsize, eps, alg::Fourier) = ceil(Int64, 0.5*(stepsize/(π*eps))^2)
terms_needed(W, stepsize, eps, alg::Wiktorsson) = begin; m=length(W); ceil(Int64, √( m*(m-1)*(m+4*(W'*W)/stepsize)/6 ) * stepsize/(2π*eps)) end

terms_needed(W, q_12::AbstractVector, stepsize, eps, alg::Fourier) = ceil(Int64, (stepsize*(q_12'*q_12)/(π*eps))^2) # here we ignored a constant factor
# In the case of the Wiktorsson generalization to Q-Wiener processes,
# we have two different non-sharp bounds. We can use the minimum of both.
function terms_needed(W, q_12::AbstractVector, stepsize, eps, alg::Wiktorsson)
    m = length(W)
    q_min, q_max = extrema(q_12)
    p1 = q_max*m*sqrt(m-1)
    p2 = sqrt(q_max*(q_12'*q_12)^3)/q_min
    return ceil(Int64, min(p1,p2)*stepsize/eps) # again ignoring a constant factor
end




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
