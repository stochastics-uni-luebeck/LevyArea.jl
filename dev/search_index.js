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
    "location": "#SRK.em",
    "page": "Home",
    "title": "SRK.em",
    "category": "function",
    "text": "em(drift::Function, diffusion::Function,\n    𝓘::AbstractRange,\n    m::Integer, X₀::AbstractVector{<:Real},\n    dW::Union{AbstractArray,Nothing}=nothing)\n\nTODO: write docs\n\n\n\n\n\n"
},

{
    "location": "#SRK.simdoubleintegrals",
    "page": "Home",
    "title": "SRK.simdoubleintegrals",
    "category": "function",
    "text": "simdoubleintegrals(W::AbstractVector{AbstractFloat64}, h::Real, C::Real=1.0)\nsimdoubleintegrals(W::Real, h::Real=1.0)\n\nSimulates an approximation of the iterated Itô-integrals int_0^hint_0^sdW_i(t)dW_j(s) for all pairs 1le i j le m of the given m-dimensional Brownian motion with step size h. For the used algorithm see [Wiktorsson2001].\n\nIn the case of a scalar Brownian motion the integral can be explicitly calculated as int_0^hint_0^sdW(t)dW(s) = frac12W(h)^2 - frac12h.\n\n\n\n\n\n"
},

{
    "location": "#SRK.sri",
    "page": "Home",
    "title": "SRK.sri",
    "category": "function",
    "text": "sri(drift::Function, diffusion::Function,\n    𝓘::AbstractRange,\n    m::Integer, X₀::AbstractVector{<:Real},\n    dW::Union{AbstractArray,Nothing}=nothing)\n\nTODO: write docs\n\n\n\n\n\n"
},

{
    "location": "#Exported-Functions-1",
    "page": "Home",
    "title": "Exported Functions",
    "category": "section",
    "text": "Modules = [SRK]\nPrivate = false\nOrder = [:type, :function]"
},

{
    "location": "#SRK.createBrownianPath-Tuple{AbstractRange,Integer}",
    "page": "Home",
    "title": "SRK.createBrownianPath",
    "category": "method",
    "text": "createBrownianPath(I::AbstractRange,m::Integer)\n\nGenerate an m-dimensional Brownian path at the timepoints specified by I. The path always starts at 0_m.\n\n\n\n\n\n"
},

{
    "location": "#SRK.restrict",
    "page": "Home",
    "title": "SRK.restrict",
    "category": "function",
    "text": "restrict(I::AbstractRange,dW::AbstractArray,k::Integer=1)\n\nRestrict a path of a Brownian motion given by a range of time points and the corresponding Brownian increments to frac12^k as many time points. The lengths of I and the first dimension of dW must match and be of the form 2^n+1 for some nge k.\n\n\n\n\n\n"
},

{
    "location": "#SRK.simdoubleintegrals_n-Tuple{AbstractArray{#s21,1} where #s21<:AbstractFloat,Integer}",
    "page": "Home",
    "title": "SRK.simdoubleintegrals_n",
    "category": "method",
    "text": "simdoubleintegrals(W::AbstractVector{AbstractFloat64}, n::Integer)\n\nSimulates an approximation of all one-time iterated Itô-integrals of the given Brownian motions with step size 1. The algorithm is taken from [Wiktorsson2001].\n\nInput:  W   the increments of m Brownian motions, where m = length(W)         n   number of terms in the approximation of the stochastic area integral Output: I[i,j] is an approximation of int_0^1W_i(s)dW_j(s)\n\n\n\n\n\n"
},

{
    "location": "#Non-Exported-Functions-1",
    "page": "Home",
    "title": "Non-Exported Functions",
    "category": "section",
    "text": "Modules = [SRK]\nPublic = false\nOrder = [:type, :function]"
},

{
    "location": "#Index-1",
    "page": "Home",
    "title": "Index",
    "category": "section",
    "text": ""
},

{
    "location": "#References-1",
    "page": "Home",
    "title": "References",
    "category": "section",
    "text": "[Wiktorsson2001]: \"Joint characteristic function and simultaneous simulation of iterated Itô integrals for multiple independent Brownian motions.\" The Annals of Applied Probability 11.2: 470-487."
},

]}
