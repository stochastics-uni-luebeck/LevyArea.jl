module LevyArea

# Imports
import LinearAlgebra: mul!
import Random: randn!, default_rng
import SpecialFunctions: trigamma


"""
    abstract type AbstractIteratedIntegralAlgorithm end

Abstract type for algorithms for the simulation of iterated integrals.

```jldoctest; setup=:(using InteractiveUtils; using LevyArea)
julia> subtypes(AbstractIteratedIntegralAlgorithm)
4-element Vector{Any}:
 Fourier
 Milstein
 MronRoe
 Wiktorsson
```
"""
abstract type AbstractIteratedIntegralAlgorithm end


# Define error criteria
include("error_norms.jl")
export MaxL2
export FrobeniusL2

# Simulate in terms of Levy Area
include("iterated_integrals.jl")
export iterated_integrals

# Properties of Levy Area algorithms
include("alg_utils.jl")
export terms_needed
export optimal_algorithm

# Levy Area algorithms
include("fourier.jl")
include("milstein.jl")
include("wiktorsson.jl")
include("mronroe.jl")
export AbstractIteratedIntegralAlgorithm
export Fourier
export Milstein
export Wiktorsson
export MronRoe

const ITER_INT_ALGS = [Fourier(),Milstein(),Wiktorsson(),MronRoe()]


end # module
