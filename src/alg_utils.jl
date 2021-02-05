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

# Default implementations in terms of the MaxL2 norm
function errcoeff(m, stepsize, alg::AbstractIteratedIntegralAlgorithm, ::FrobeniusL2)
    return √(m^2-m) * errcoeff(m, stepsize, alg, MaxL2())
end
function errcoeff(m, q_12, stepsize, alg::AbstractIteratedIntegralAlgorithm, ::MaxL2)
    maxqq = maximum(q_12[i]*q_12[j] for i=1:m for j=1:i-1)
    return maxqq * errcoeff(m, stepsize, alg, MaxL2())
end
function errcoeff(m, q_12, stepsize, alg::AbstractIteratedIntegralAlgorithm, ::FrobeniusL2)
    trQ_sq = abs2(sum(abs2, q_12))
    tr_Qsq = sum(abs2∘abs2, q_12)
    return √(trQ_sq-tr_Qsq) * errcoeff(m, stepsize, alg, MaxL2())
end

"""
    norv(dim, n, alg::AbstractIteratedIntegralAlgorithm)

Returns the number of random numbers needed to simulate the iterated integrals
for a Wiener process of dimension `dim` and with truncation parameter `n`.
"""
norv(dim, n, alg::AbstractIteratedIntegralAlgorithm)

"""
    effective_cost(dim, stepsize, eps, alg, norm)

Returns the number of random numbers needed to simulate the iterated integrals
with the given parameters to the given precision.
"""
function effective_cost(dim, stepsize, eps, alg, norm)
    norv(m, terms_needed(dim, stepsize, eps, alg, norm), alg)
end