"""
    convorder(alg::AbstractIteratedIntegralAlgorithm)

Returns the convergence order of the algorithm w.r.t. the truncation parameter.

See also: [`errcoeff`](@ref)
"""
function convorder end

"""
    errcoeff(dim, stepsize, alg, norm)
    errcoeff(dim, q_12, stepsize, alg, norm)

Returns the coefficient of the truncation parameter in the error estimate.
I.e. the error estimate is of the form

```math
\\lVert I(h)-\\tilde{I}^{(p)}(h) \\rVert_* ≤ \\mathrm{errcoeff}(m,h) \\cdot p^{-γ}
```

where the norm is given by `norm`, the approximation ``\\tilde{I}^{(p)}`` is 
calculated using `alg` and ``γ`` is the order of convergence given by [`convorder(alg)`](@ref).
"""
function errcoeff end

"""
    terms_needed(dim, stepsize, eps, alg, norm)
    
Returns the number of terms in the approximating sum that is needed to ensure an error
in the given norm of at most `eps`.
This depends on the dimension of the Wiener process `dim`, the current stepsize and the chosen algorithm.
    
See also: [`AbstractIteratedIntegralAlgorithm`](@ref), [`AbstractErrorNorm`](@ref)

# Examples
```jldoctest; setup=:(using LevyArea)
julia> h = 1/128;

julia> terms_needed(10, h, h^(3/2), Milstein(), MaxL2())
7
```

# Implementation
New algorithms should only have to implement [`errcoeff`](@ref) and [`convorder`](@ref).
"""
function terms_needed(dim, stepsize, eps, alg::AbstractIteratedIntegralAlgorithm, norm::AbstractErrorNorm)
    ceil(Int64, (errcoeff(dim, stepsize, alg, norm)/eps)^(1//convorder(alg)) )
end

"""
    terms_needed(dim, q_12, stepsize, eps, alg, norm)
    
Used for finite-dimensional approximations of a Q-Wiener process with covariance matrix
``Q = Q^\\frac{1}{2}*Q^\\frac{1}{2}``. Here `q_12` is a vector of the eigenvalues of ``Q^\\frac{1}{2}``;
the square root of the covariance matrix. Equivalently these are the square roots of the eigenvalues of ``Q``.

# Examples
```jldoctest; setup=:(using LevyArea)
julia> h = 1/128;

julia> dim = 10;

julia> q = [1/k^2 for k=1:dim];

julia> terms_needed(dim, sqrt.(q), h, h^(3/2), Milstein(), FrobeniusL2())
9
```
"""
function terms_needed(dim, q_12, stepsize, eps, alg::AbstractIteratedIntegralAlgorithm, norm::AbstractErrorNorm)
    length(q_12) == dim || throw(ArgumentError("Length of q_12 must be equal to the dimension."))
    ceil(Int64, (errcoeff(dim, q_12, stepsize, alg, norm)/eps)^(1//convorder(alg)) )
end

"""
    norv(dim, n, alg::AbstractIteratedIntegralAlgorithm)

Returns the number of random numbers needed to simulate the iterated integrals
for a Wiener process of dimension `dim` and with truncation parameter `n`.
"""
norv(dim, n, alg::AbstractIteratedIntegralAlgorithm)

"""
    effective_cost(dim, stepsize, eps, alg, norm)
    effective_cost(dim, q_12, stepsize, eps, alg, norm)

Returns the number of random numbers needed to simulate the iterated integrals
with the given parameters.
"""
function effective_cost(dim, stepsize, eps, alg, norm)
    norv(dim, terms_needed(dim, stepsize, eps, alg, norm), alg)
end
function effective_cost(dim, q_12, stepsize, eps, alg, norm)
    norv(dim, terms_needed(dim, q_12, stepsize, eps, alg, norm), alg)
end


"""
    optimal_algorithm(dim, stepsize, eps=stepsize^(3/2), norm=MaxL2())
    optimal_algorithm(dim, q_12, stepsize, eps, norm=FrobeniusL2())

Returns the optimal algorithm for the given parameters,
i.e. the algorithm that needs to simulate the fewest random numbers
to achieve the desired precision `eps`.

# Examples
```jldoctest; setup=:(using LevyArea)
julia> h = 1/128;

julia> optimal_algorithm(10, h, h^(3/2), MaxL2())
MR()

julia> optimal_algorithm(10, 1.0./(1:10).^2, h, h^(3/2), FrobeniusL2())
Milstein()
```
"""
function optimal_algorithm(dim, stepsize, eps=stepsize^(3/2), norm::AbstractErrorNorm=MaxL2())
    ind = argmin([effective_cost(dim, stepsize, eps, alg, norm) for alg ∈ ITER_INT_ALGS])
    return ITER_INT_ALGS[ind]
end
function optimal_algorithm(dim, q_12, stepsize, eps, norm::AbstractErrorNorm=FrobeniusL2())
    ind = argmin([effective_cost(dim, q_12, stepsize, eps, alg, norm) for alg ∈ ITER_INT_ALGS])
    return ITER_INT_ALGS[ind]
end
