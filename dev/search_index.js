var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#SRK.jl-1",
    "page": "Home",
    "title": "SRK.jl",
    "category": "section",
    "text": "Stochastic Runge-Kutta methods in JuliaThis package tries to implement state-of-the-art SRK methods for solving stochastic differential equations (SDEs)."
},

{
    "location": "#SRK.simdoubleintegrals-Tuple{Array{Float64,1},Integer}",
    "page": "Home",
    "title": "SRK.simdoubleintegrals",
    "category": "method",
    "text": "simdoubleintegrals(W::Vector{Float64}, n::Integer)\n\nSimulates an approximation of all one-time iterated Itô-integrals of the given Brownian motions with step size 1. The algorithm is taken from [Wiktorsson2001].\n\nInput:  W   the increments of m Brownian motions, where m = length(W)         n   number of terms in the approximation of the stochastic area integral Output: I[i,j] is an approximation of int_0^hW_i(s)dW_j(s)\n\nReferences\n\n[Wiktorsson2001]: \"Joint characteristic function and simultaneous simulation of iterated Itô integrals for multiple independent Brownian motions.\" The Annals of Applied Probability 11.2: 470-487.\n\n\n\n\n\n"
},

{
    "location": "#Functions-1",
    "page": "Home",
    "title": "Functions",
    "category": "section",
    "text": "Modules = [SRK]\nPrivate = false\nOrder = [:type, :function]"
},

{
    "location": "#Index-1",
    "page": "Home",
    "title": "Index",
    "category": "section",
    "text": ""
},

]}
