module SRK

import LinearAlgebra: mul!, lmul!, I
import Random: GLOBAL_RNG, randn!
import SpecialFunctions: trigamma

# remove (only needed for 'old' algorithm)
#import LinearAlgebra: kron, Diagonal, tril!
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
