module SRK

import LinearAlgebra: kron, lmul!, I, Diagonal
import SparseArrays: sparse
import SpecialFunctions: trigamma


export simdoubleintegrals


include("simdoubleintegrals.jl")


end # module
