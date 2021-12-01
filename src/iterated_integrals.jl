"""
    ito_correction!(I, h=1)

Applies the Itô-correction for iterated integrals to `I`.
This amounts to subtracting ``\\frac{1}{2}h`` from every element on the diagonal.

# Example
```jldoctest; setup=:(using LevyArea)
julia> M = ones(5,5);

julia> LevyArea.ito_correction!(M, 0.5)


julia> M
5×5 Matrix{Float64}:
 0.75  1.0   1.0   1.0   1.0
 1.0   0.75  1.0   1.0   1.0
 1.0   1.0   0.75  1.0   1.0
 1.0   1.0   1.0   0.75  1.0
 1.0   1.0   1.0   1.0   0.75
```

"""
function ito_correction!(I, h=1)
    m,n = size(I)
    m == n || throw(DimensionMismatch("Matrix is not square: dimensions are $(size(I))"))
    @inbounds for i=1:m
        I[i,i] -= h/2
    end
end

"""
    iterated_integrals(W::AbstractVector, h, eps=h^(3/2);
        ito_correction=true,
        error_norm=MaxL2(),
        alg=optimal_algorithm(length(W),h,eps,error_norm)
    )

Simulates an approximation of the iterated stochastic integrals
``\\int_0^h\\int_0^sdW_i(t)dW_j(s)`` for all pairs ``1\\le i, j \\le m``
of the given ``m``-dimensional Brownian motion with step size h.

# Examples
```jldoctest; setup=:(using LinearAlgebra; using LevyArea)
julia> h = 1/2;

julia> W = [1.0, 0.5]
2-element Vector{Float64}:
 1.0
 0.5

julia> diag(iterated_integrals(W, h, h^(3/2))) ≈ 0.5*W.^2 .- 0.5h
true
```
"""
function iterated_integrals(W::AbstractVector{T}, h::Real, eps::Real=h^(3/2);
    ito_correction=true,
    error_norm::AbstractErrorNorm=MaxL2(),
    alg::AbstractIteratedIntegralAlgorithm=optimal_algorithm(length(W),h,eps,error_norm)
) where {T<:AbstractFloat}
    m = length(W)
    n = terms_needed(m, h, eps, alg, error_norm)
    I = levyarea(W/√h, n, alg)
    if ito_correction
        ito_correction!(I)
    end
    I .= 0.5.*W.*W' .+ h.*I
    return I
end

"""
    iterated_integrals(W::AbstractVector, q_12::AbstractVector, h, eps; 
        ito_correction=true,
        error_norm=FrobeniusL2(),
        alg=optimal_algorithm(length(W),q_12,h,eps,error_norm)
    )

Simulates an approximation of the iterated stochastic integrals for finite-dimensional approximations of
a Q-Wiener process with covariance matrix ``Q = Q^\\frac{1}{2}*Q^\\frac{1}{2}``.
Here `q_12` is a vector of the eigenvalues of ``Q^\\frac{1}{2}``; the square root of the covariance matrix.
Equivalently these are the square roots of the eigenvalues of ``Q``.

# Examples
```jldoctest; setup=:(using LinearAlgebra; using LevyArea)
julia> h = 0.01; dim=10; q = [1/k^2 for k=1:dim];

julia> W = √h * sqrt.(q) .* randn(dim);

julia> diag(iterated_integrals(W,sqrt.(q),h,h^(3/2))) ≈ 0.5*W.^2 .- 0.5*h*q
true
```
"""
function iterated_integrals(W::AbstractVector{T}, q_12::AbstractVector, h::Real, eps::Real;
    ito_correction=true,
    error_norm::AbstractErrorNorm=FrobeniusL2(),
    alg::AbstractIteratedIntegralAlgorithm=optimal_algorithm(length(W),q_12,h,eps,error_norm)
) where {T<:AbstractFloat}
    m = length(W)
    n = terms_needed(m, q_12, h, eps, alg, error_norm)
    I = levyarea(W./q_12./√h, n, alg)
    if ito_correction
        ito_correction!(I)
    end
    I .= 0.5.*W.*W' .+ h.*q_12'.*I.*q_12 # scale correctly
    return I
end

"""
    iterated_integrals(W::Real, h::Real, eps::Real=0.0; ito_correction=true, kwargs...)

In the case of a scalar Brownian motion the integral can be explicitly
calculated as ``\\int_0^h\\int_0^sdW(t)dW(s) = \\frac{1}{2}W(h)^2 - \\frac{1}{2}h``.

The parameter `eps` (as well as all additional keyword arguments) has no effect but is available 
to provide the same interface as the multidimensional version.
"""
iterated_integrals(W::Real, h::Real, eps::Real=0.0; ito_correction=true, kwargs...) = ito_correction ? 0.5W^2 - 0.5h : 0.5W^2

"""
    iterated_integrals(W::Real, q_12::Real, h::Real, eps::Real; ito_correction=true, kwargs...)

In the case of a scalar Q-Wiener process with (scalar) covariance Q the integral can be explicitly
calculated as ``\\int_0^h\\int_0^sdW(t)dW(s) = \\frac{1}{2}W(h)^2 - \\frac{1}{2}hQ``.

Note that, as in the multidimensional case, the parameter `q_12` denotes the square root of the covariance.

The parameter `eps` (as well as all additional keyword arguments) has no effect but is available 
to provide the same interface as the multidimensional version.
"""
iterated_integrals(W::Real, q_12::Real, h::Real, eps::Real; ito_correction=true, kwargs...) = ito_correction ? 0.5W^2 - 0.5*h*q_12^2 : 0.5W^2