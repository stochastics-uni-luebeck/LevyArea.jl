module IteratedIntegrals

# Imports
import LinearAlgebra: mul!, lmul!, BLAS
import Random: GLOBAL_RNG, randn!
import SpecialFunctions: trigamma

# remove (only needed for 'old' algorithm)
#import LinearAlgebra: kron, Diagonal, tril!, I
#import SparseArrays: sparse



# Exports
export simdoubleintegrals
export sri
export em


abstract type AbstractIteratedIntegralAlgorithm end
struct Milstein <: AbstractIteratedIntegralAlgorithm end
struct Wiktorsson <: AbstractIteratedIntegralAlgorithm end



include("simdoubleintegrals.jl")
include("sri.jl")
include("em.jl")
include("utils.jl")



end # module
