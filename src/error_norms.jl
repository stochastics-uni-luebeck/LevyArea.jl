# Note: All the estimates and inequalities below only apply
# for our iterated integrals because the Levy area and the 
# corresponding error terms are skew-symmetric!

"""
    abstract type AbstractErrorNorm end

Abstract type for different kind of errors one might consider.
The most important are `MaxL2` and `FrobeniusL2`. These are actually
special cases of the `MaxLp{p}` and `SchattenqLp{p,q}` norms:
```julia
const MaxL2 = MaxLp{2}
const FrobeniusL2 = SchattenqLp{2,2}
```

All currently defined norms:
```jldoctest; setup=:(using InteractiveUtils; using IteratedIntegrals)
julia> subtypes(IteratedIntegrals.AbstractErrorNorm)
2-element Vector{Any}:
 IteratedIntegrals.LpMax
 IteratedIntegrals.MaxLp
 IteratedIntegrals.SchattenqLp
```
"""
abstract type AbstractErrorNorm end


# maximum of entry-wise Lp norms
# this requires new Lp convergence proofs
# and should be defined for each new algorithm
struct MaxLp{p} <: AbstractErrorNorm end

function errcoeff(m, q_12, stepsize, alg::AbstractIteratedIntegralAlgorithm, ::MaxLp{p}) where {p}
    maxqq = maximum(q_12[i]*q_12[j] for i=1:m for j=1:i-1)
    return maxqq * errcoeff(m, stepsize, alg, MaxLp{p}())
end

# lq norm of eigenvalues of matrix of entry-wise Lp norms
# (= Schatten-q norm of matrix of entry-wise Lp norms)
# we can get this from the MaxLp norms
# p=q=2 is the FrobeniusL2 norm
struct SchattenqLp{p,q} <: AbstractErrorNorm end

function errcoeff(m, stepsize, alg::AbstractIteratedIntegralAlgorithm, ::SchattenqLp{p,q}) where {p,q}
    return (m^2-m)^(1/q) * errcoeff(m, stepsize, alg, MaxLp{p}())
end
function errcoeff(m, q_12, stepsize, alg::AbstractIteratedIntegralAlgorithm, ::SchattenqLp{p,q}) where {p,q}
    trQq_sq = abs2(sum(x->x^q, q_12))
    tr_Qqsq = sum(x->x^2q, q_12)
    return (trQq_sq-tr_Qqsq)^(1/q) * errcoeff(m, stepsize, alg, MaxLp{p}())
end

# Lp norm of maximum of matrix
# can be bounded by the SchattenqLp norm with q=p
# (this relies on the fact that Schatten-p-Lp == Lp-Schatten-p)
# the factor 1/2^1/p comes from the antisymmetry
struct LpMax{p} <: AbstractErrorNorm end

function errcoeff(m, stepsize, alg::AbstractIteratedIntegralAlgorithm, ::LpMax{p}) where {p}
    return 1/(2^(1/p)) * errcoeff(m, stepsize, alg, SchattenqLp{p,p}())
end
function errcoeff(m, q_12, stepsize, alg::AbstractIteratedIntegralAlgorithm, ::LpMax{p}) where {p}
    return 1/(2^(1/p)) * errcoeff(m, q_12, stepsize, alg, SchattenqLp{p,p}())
end


# most important/common norms
const MaxL2 = MaxLp{2}
const FrobeniusL2 = SchattenqLp{2,2}