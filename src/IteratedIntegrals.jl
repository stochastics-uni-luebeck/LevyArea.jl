module IteratedIntegrals

# Imports
import LinearAlgebra: mul!, lmul!, BLAS
import Random: randn!
import SpecialFunctions: trigamma

# remove (only needed for 'old' algorithm)
#import LinearAlgebra: kron, Diagonal, tril!, I
#import SparseArrays: sparse
#import Random: GLOBAL_RNG

# Types
abstract type AbstractIteratedIntegralAlgorithm end
struct Fourier <: AbstractIteratedIntegralAlgorithm end
struct Wiktorsson <: AbstractIteratedIntegralAlgorithm end

# Exports
export simiterintegrals
export sri
export em

export AbstractIteratedIntegralAlgorithm
export Fourier
export Wiktorsson




include("simiterintegrals.jl")
include("sri.jl")
include("em.jl")
include("utils.jl")



end # module
