# LevyArea.jl
*Iterated Stochastic Integrals in Julia*

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://stochastics-uni-luebeck.github.io/LevyArea.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://stochastics-uni-luebeck.github.io/LevyArea.jl/dev)
[![Build Status](https://github.com/stochastics-uni-luebeck/LevyArea.jl/workflows/CI/badge.svg)](https://github.com/stochastics-uni-luebeck/LevyArea.jl/actions)

This package implements state-of-the-art methods for the simulation of iterated stochastic integrals.
These appear e.g. in higher order algorithms for the solution of stochastic (partial) differential equations.

## Installation

This package can be installed from the Julia package manager (type <kbd>]</kbd>)
```julia
pkg> add LevyArea
```

## Usage Example

```julia
julia> using LevyArea

julia> h = 1/100;

julia> W = sqrt(h) * randn(5)
5-element Array{Float64,1}:
 -0.15140313307420128
 -0.031565386565872114
  0.04288593819444138
 -0.03478687621740065
  0.07116134579533281

julia> iterated_integrals(W, h, h^(3/2))
5×5 Array{Float64,2}:
  0.00646145  -0.00288169  -0.00092251    0.00208285  -0.00912994
  0.00766079  -0.00450181  -0.00126378    0.00517929  -0.00967469
 -0.00557056  -8.99281e-5  -0.0040804    -0.00165624   0.00434686
  0.003184    -0.00408123   0.000164376  -0.00439494  -0.00381448
 -0.00164411   0.00742845  -0.00129504    0.001339    -0.00246803

julia> 0.5*W.^2 .- 0.5h
5-element Array{Float64,1}:
  0.006461454352342151
 -0.00450181318547353
 -0.004080398152591277
 -0.004394936621517622
 -0.0024680314322985345
```
