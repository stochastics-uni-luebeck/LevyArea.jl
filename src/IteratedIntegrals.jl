module IteratedIntegrals

# Imports
import LinearAlgebra: mul!, lmul!, BLAS
import Random: randn!, default_rng
import SpecialFunctions: trigamma

# Types
"""
    abstract type AbstractIteratedIntegralAlgorithm end

Abstract type for algorithms for the simulation of iterated integrals.

```jldoctest; setup=:(using InteractiveUtils; using IteratedIntegrals)
julia> subtypes(AbstractIteratedIntegralAlgorithm)
4-element Vector{Any}:
 Fourier
 MR
 Milstein
 Wiktorsson
```
"""
abstract type AbstractIteratedIntegralAlgorithm end
"""
    abstract type AbstractErrorNorm end

Abstract type for different kind of errors one might consider.

```jldoctest; setup=:(using InteractiveUtils; using IteratedIntegrals)
julia> subtypes(IteratedIntegrals.AbstractErrorNorm)
2-element Vector{Any}:
 FrobeniusL2
 MaxL2
```
"""
abstract type AbstractErrorNorm end
struct MaxL2 <: AbstractErrorNorm end
struct FrobeniusL2 <: AbstractErrorNorm end

# Exports
export simiterintegrals
export terms_needed
export optimal_algorithm

export MaxL2
export FrobeniusL2

export AbstractIteratedIntegralAlgorithm
export Fourier
export Milstein
export Wiktorsson
export MR


# Simulate in terms of Levy Area
include("simiterintegrals.jl")
# Properties of Levy Area algorithms
include("alg_utils.jl")
# Levy Area algorithms
include("fourier.jl")
include("milstein.jl")
include("wiktorsson.jl")
include("mr.jl")

const ITER_INT_ALGS = [Fourier(),Milstein(),Wiktorsson(),MR()]

# Other stuff
include("sri.jl")
include("em.jl")
include("utils.jl")



end # module
