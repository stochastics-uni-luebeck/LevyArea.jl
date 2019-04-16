module IteratedIntegrals

import LinearAlgebra: mul!, lmul!, BLAS
import Random: GLOBAL_RNG, randn!
import SpecialFunctions: trigamma

# remove (only needed for 'old' algorithm)
#import LinearAlgebra: kron, Diagonal, tril!, I
#import SparseArrays: sparse
#

export simdoubleintegrals
export sri
export em


include("simdoubleintegrals.jl")
include("sri.jl")
include("em.jl")
include("utils.jl")



end # module
