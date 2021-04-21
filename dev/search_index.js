var documenterSearchIndex = {"docs":
[{"location":"functions/#Function-reference","page":"Index","title":"Function reference","text":"","category":"section"},{"location":"functions/","page":"Index","title":"Index","text":"DocTestSetup = :(using LinearAlgebra; using IteratedIntegrals)","category":"page"},{"location":"functions/#Main-Functions","page":"Index","title":"Main Functions","text":"","category":"section"},{"location":"functions/","page":"Index","title":"Index","text":"simiterintegrals\nIteratedIntegrals.levyarea\nIteratedIntegrals.AbstractIteratedIntegralAlgorithm\nIteratedIntegrals.AbstractErrorNorm\noptimal_algorithm","category":"page"},{"location":"functions/#IteratedIntegrals.simiterintegrals","page":"Index","title":"IteratedIntegrals.simiterintegrals","text":"simiterintegrals(W::AbstractVector, h, eps=h^(3/2); ito_correction=true, alg=MR())\n\nSimulates an approximation of the iterated stochastic integrals int_0^hint_0^sdW_i(t)dW_j(s) for all pairs 1le i j le m of the given m-dimensional Brownian motion with step size h.\n\nExamples\n\njulia> h = 1/2;\n\njulia> W = [1.0, 0.5]\n2-element Vector{Float64}:\n 1.0\n 0.5\n\njulia> diag(simiterintegrals(W, h, h^(3/2))) ≈ 0.5*W.^2 .- 0.5h\ntrue\n\n\n\n\n\nsimiterintegrals(W::AbstractVector, q_12::AbstractVector, h, eps; ito_correction=true, alg=MR())\n\nSimulates an approximation of the iterated stochastic integrals for finite-dimensional approximations of a Q-Wiener process with covariance matrix Q = Q^frac12*Q^frac12. Here q_12 is a vector of the eigenvalues of Q^frac12; the square root of the covariance matrix. Equivalently these are the square roots of the eigenvalues of Q.\n\n\n\n\n\nsimiterintegrals(W::Real, h::Real=1.0, eps::Real=0.0; ito_correction=true, kwargs...)\n\nIn the case of a scalar Brownian motion the integral can be explicitly calculated as int_0^hint_0^sdW(t)dW(s) = frac12W(h)^2 - frac12h.\n\nThe parameter eps (as well as all additional keyword arguments) has no effect but is available  to provide the same interface as the multidimensional version.\n\n\n\n\n\n","category":"function"},{"location":"functions/#IteratedIntegrals.levyarea","page":"Index","title":"IteratedIntegrals.levyarea","text":"levyarea(W, n, alg::Fourier)\n\nSimulates an approximation of the iterated Itô-integrals int_0^1W_sotimes dW_s of the given m-dimensional increment of a Wiener process with step size 1. The parameter n specifies the number of terms in the approximation and thus determines the accuracy. This algorithm is based on a Fourier expansion of the Wiener process. The algorithm needs approximately 2cdot m^2+2cdot mcdot n Float's and 2cdot mcdot n random numbers. The time complexity is mathcalO(m^2cdot n).\n\n\n\n\n\nlevyarea(W, n, alg::Milstein)\n\nSimulates an approximation of the iterated Itô-integrals int_0^1W_sotimes dW_s of the given m-dimensional increment of a Wiener process with step size 1. The parameter n specifies the number of terms in the approximation and thus determines the accuracy. This is an efficient implementation of the algorithm proposed in Milstein, 1994. It is based on a Fourier expansion of the Wiener process. The algorithm needs approximately 2cdot m^2+2cdot mcdot n+m Float64's. The time complexity is mathcalO(m^2cdot n).\n\n\n\n\n\nlevyarea(W, n, alg::Wiktorsson)\n\nSimulates an approximation of the iterated Itô-integrals int_0^1W_sotimes dW_s of the given m-dimensional increment of a Wiener process with step size 1. The parameter n specifies the number of terms in the approximation and thus determines the accuracy. This is an efficient implementation of the algorithm proposed in Wiktorsson, 2001. It is based on the Fourier method from Milstein but incorporates an additional tail sum approximation. The algorithm needs approximately 2cdot m^2+2cdot mcdot n+m Float64's. The time complexity is mathcalO(m^2cdot n).\n\n\n\n\n\nlevyarea(W, n, alg::MR)\n\nSimulates an approximation of the iterated Itô-integrals int_0^1W_sotimes dW_s of the given m-dimensional increment of a Wiener process with step size 1. The parameter n specifies the number of terms in the approximation and thus determines the accuracy. This is an efficient implementation of the algorithm proposed in Mrongowius & Rößler, 2021. It is based on the Fourier method from Milstein but incorporates an improved tail sum approximation. The algorithm needs approximately m^2+2cdot mcdot n Float's  and 12m^2+2cdot mcdot n + 12m random numbers. The time complexity is mathcalO(m^2cdot n).\n\n\n\n\n\n","category":"function"},{"location":"functions/#IteratedIntegrals.AbstractIteratedIntegralAlgorithm","page":"Index","title":"IteratedIntegrals.AbstractIteratedIntegralAlgorithm","text":"abstract type AbstractIteratedIntegralAlgorithm end\n\nAbstract type for algorithms for the simulation of iterated integrals.\n\njulia> subtypes(AbstractIteratedIntegralAlgorithm)\n4-element Vector{Any}:\n Fourier\n MR\n Milstein\n Wiktorsson\n\n\n\n\n\n","category":"type"},{"location":"functions/#IteratedIntegrals.AbstractErrorNorm","page":"Index","title":"IteratedIntegrals.AbstractErrorNorm","text":"abstract type AbstractErrorNorm end\n\nAbstract type for different kind of errors one might consider.\n\njulia> subtypes(IteratedIntegrals.AbstractErrorNorm)\n2-element Vector{Any}:\n FrobeniusL2\n MaxL2\n\n\n\n\n\n","category":"type"},{"location":"functions/#IteratedIntegrals.optimal_algorithm","page":"Index","title":"IteratedIntegrals.optimal_algorithm","text":"optimal_algorithm(dim, stepsize, eps=stepsize^(3/2), norm=MaxL2())\noptimal_algorithm(dim, q_12, stepsize, eps, norm=FrobeniusL2())\n\nReturns the optimal algorithm for the given parameters, i.e. the algorithm that needs to simulate the fewest random numbers.\n\nExamples\n\njulia> h = 1/128;\n\njulia> optimal_algorithm(10, h, h^(3/2), MaxL2())\nMR()\n\njulia> optimal_algorithm(10, 1.0./(1:10).^2, h, h^(3/2), FrobeniusL2())\nMilstein()\n\n\n\n\n\n","category":"function"},{"location":"functions/#Algorithmic-properties","page":"Index","title":"Algorithmic properties","text":"","category":"section"},{"location":"functions/","page":"Index","title":"Index","text":"terms_needed\nIteratedIntegrals.convorder\nIteratedIntegrals.errcoeff\nIteratedIntegrals.norv\nIteratedIntegrals.effective_cost","category":"page"},{"location":"functions/#IteratedIntegrals.terms_needed","page":"Index","title":"IteratedIntegrals.terms_needed","text":"terms_needed(dim, stepsize, eps, alg, norm)\n\nReturns the number of terms in the approximating sum that is needed to ensure an error in the given norm of at most eps. This depends on the dimension of the Wiener process dim, the current stepsize and the chosen algorithm.\n\nSee also: AbstractIteratedIntegralAlgorithm, AbstractErrorNorm\n\nExamples\n\njulia> h = 1/128;\n\njulia> terms_needed(10, h, h^(3/2), Milstein(), MaxL2())\n7\n\nImplementation\n\nNew algorithms should only have to implement errcoeff and convorder.\n\n\n\n\n\nterms_needed(dim, q_12, stepsize, eps, alg, norm)\n\nUsed for finite-dimensional approximations of a Q-Wiener process with covariance matrix Q = Q^frac12*Q^frac12. Here q_12 is a vector of the eigenvalues of Q^frac12; the square root of the covariance matrix. Equivalently these are the square roots of the eigenvalues of Q.\n\nExamples\n\njulia> h = 1/128;\n\njulia> dim = 10;\n\njulia> q = [1/k^2 for k=1:dim];\n\njulia> terms_needed(dim, sqrt.(q), h, h^(3/2), Milstein(), FrobeniusL2())\n9\n\n\n\n\n\n","category":"function"},{"location":"functions/#IteratedIntegrals.convorder","page":"Index","title":"IteratedIntegrals.convorder","text":"convorder(alg::AbstractIteratedIntegralAlgorithm)\n\nReturns the convergence order of the algorithm w.r.t. the truncation parameter.\n\nSee also: errcoeff\n\n\n\n\n\n","category":"function"},{"location":"functions/#IteratedIntegrals.errcoeff","page":"Index","title":"IteratedIntegrals.errcoeff","text":"errcoeff(dim, stepsize, alg, norm)\nerrcoeff(dim, q_12, stepsize, alg, norm)\n\nReturns the coefficient of the truncation parameter in the error estimate. I.e. the error estimate is of the form\n\nlVert I(h)-tildeI^(p)(h) rVert_*  mathrmerrcoeff(mh) cdot p^-γ\n\nwhere the norm is given by norm, the approximation tildeI^(p) is  calculated using alg and γ is the order of convergence given by convorder(alg).\n\n\n\n\n\n","category":"function"},{"location":"functions/#IteratedIntegrals.norv","page":"Index","title":"IteratedIntegrals.norv","text":"norv(dim, n, alg::AbstractIteratedIntegralAlgorithm)\n\nReturns the number of random numbers needed to simulate the iterated integrals for a Wiener process of dimension dim and with truncation parameter n.\n\n\n\n\n\n","category":"function"},{"location":"functions/#IteratedIntegrals.effective_cost","page":"Index","title":"IteratedIntegrals.effective_cost","text":"effective_cost(dim, stepsize, eps, alg, norm)\neffective_cost(dim, q_12, stepsize, eps, alg, norm)\n\nReturns the number of random numbers needed to simulate the iterated integrals with the given parameters.\n\n\n\n\n\n","category":"function"},{"location":"functions/#Helper-functions","page":"Index","title":"Helper functions","text":"","category":"section"},{"location":"functions/","page":"Index","title":"Index","text":"IteratedIntegrals.ito_correction!","category":"page"},{"location":"functions/#IteratedIntegrals.ito_correction!","page":"Index","title":"IteratedIntegrals.ito_correction!","text":"ito_correction!(I, h=1)\n\nApplies the Itô-correction for iterated integrals to I. This amounts to subtracting frac12h from every element on the diagonal.\n\nExample\n\njulia> M = ones(5,5);\n\njulia> IteratedIntegrals.ito_correction!(M, 0.5)\n\n\njulia> M\n5×5 Matrix{Float64}:\n 0.75  1.0   1.0   1.0   1.0\n 1.0   0.75  1.0   1.0   1.0\n 1.0   1.0   0.75  1.0   1.0\n 1.0   1.0   1.0   0.75  1.0\n 1.0   1.0   1.0   1.0   0.75\n\n\n\n\n\n","category":"function"},{"location":"functions/#Index","page":"Index","title":"Index","text":"","category":"section"},{"location":"functions/","page":"Index","title":"Index","text":"","category":"page"},{"location":"references/#References","page":"References","title":"References","text":"","category":"section"},{"location":"references/#milstein1994","page":"References","title":"Milstein, 1994","text":"","category":"section"},{"location":"references/","page":"References","title":"References","text":"Milstein, G. N., \"Numerical integration of stochastic differential equations.\"Vol. 313. Springer Science & Business Media, 1994.","category":"page"},{"location":"references/#mr2021","page":"References","title":"Mrongowius & Rößler, 2021","text":"","category":"section"},{"location":"references/","page":"References","title":"References","text":"Mrongowius, J. and Rößler, A., \"On the Approximation and Simulation of Iterated Stochastic Integrals and the Corresponding Lévy Areas in Terms of a Multidimensional Brownian Motion.\"Preprint, arXiv:2101.09542","category":"page"},{"location":"references/#wiktorsson2001","page":"References","title":"Wiktorsson, 2001","text":"","category":"section"},{"location":"references/","page":"References","title":"References","text":"Wiktorsson, M., \"Joint Characteristic Function and Simultaneous Simulation of Iterated Itô Integrals for Multiple Independent Brownian Motions.\"The Annals of Applied Probability 11.2: 470-487, 2001.","category":"page"},{"location":"#IteratedIntegrals.jl","page":"Home","title":"IteratedIntegrals.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Iterated Stochastic Integrals in Julia","category":"page"},{"location":"","page":"Home","title":"Home","text":"This package implements state-of-the-art methods for the simulation of iterated stochastic integrals. These appear e.g. in higher order algorithms for the solution of stochastic (partial) differential equations.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Since this package isn't registered yet, you have to use the GitHub URL of the repository:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> ]add https://github.com/fkastner/IteratedIntegrals.jl","category":"page"},{"location":"#Usage-Example","page":"Home","title":"Usage Example","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"julia> using IteratedIntegrals\n\njulia> h = 1/100;\n\njulia> W = sqrt(h) * randn(5)\n5-element Array{Float64,1}:\n -0.15140313307420128\n -0.031565386565872114\n  0.04288593819444138\n -0.03478687621740065\n  0.07116134579533281\n\njulia> simiterintegrals(W,h)\n5×5 Array{Float64,2}:\n  0.00646145  -0.00288169  -0.00092251    0.00208285  -0.00912994\n  0.00766079  -0.00450181  -0.00126378    0.00517929  -0.00967469\n -0.00557056  -8.99281e-5  -0.0040804    -0.00165624   0.00434686\n  0.003184    -0.00408123   0.000164376  -0.00439494  -0.00381448\n -0.00164411   0.00742845  -0.00129504    0.001339    -0.00246803\n\njulia> 0.5*W.^2 .- 0.5h\n5-element Array{Float64,1}:\n  0.006461454352342151\n -0.00450181318547353\n -0.004080398152591277\n -0.004394936621517622\n -0.0024680314322985345","category":"page"}]
}
