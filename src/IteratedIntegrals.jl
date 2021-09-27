module IteratedIntegrals

# Imports
import LinearAlgebra: mul!
import Random: randn!, default_rng
import SpecialFunctions: trigamma


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


# Define error criteria
include("error_norms.jl")
export MaxL2
export FrobeniusL2

# Simulate in terms of Levy Area
include("simiterintegrals.jl")
export simiterintegrals

# Properties of Levy Area algorithms
include("alg_utils.jl")
export terms_needed
export optimal_algorithm

# Levy Area algorithms
include("fourier.jl")
include("milstein.jl")
include("wiktorsson.jl")
include("mr.jl")
export AbstractIteratedIntegralAlgorithm
export Fourier
export Milstein
export Wiktorsson
export MR

const ITER_INT_ALGS = [Fourier(),Milstein(),Wiktorsson(),MR()]


end # module
