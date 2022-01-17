# LevyArea.jl
*Iterated Stochastic Integrals in Julia*

This package implements state-of-the-art methods for the simulation of iterated stochastic integrals.
These appear e.g. in higher order algorithms for the solution of stochastic (partial) differential equations.

## Installation

This package can be installed from the Julia package manager (type `]`)
```julia-repl
pkg> add LevyArea
```

## Usage Example

Load the package and generate a Wiener increment:
```julia-repl
julia> using LevyArea
julia> m = 5; # dimension of Wiener process
julia> h = 0.01; # step size or length of time interval
julia> err = 0.05; # error bound
julia> W = sqrt(h) * randn(m); # increment of Wiener process
```
Here, $W$ is the $m$-dimensional vector of increments of the driving
Wiener process on some time interval of length $h$.

The default call uses h^(3/2) as the precision and chooses the best algorithm automatically:
```julia-repl
julia> II = iterated_integrals(W,h)
```
If not stated otherwise, the default error criterion is the $\max,L^2$-error
and the function returns the $m \times m$ matrix `II` containing a realisation
of the approximate iterated stochastic integrals that correspond
to the given increment $W$.

The desired precision can be optionally provided
using a third positional argument:
```julia-repl
julia> II = iterated_integrals(W,h,err)
```
Again, the software package automatically chooses the optimal
algorithm.

To determine which algorithm is chosen by the package without simulating any iterated
stochastic integrals yet, the function `optimal_algorithm` can
be used. The arguments to this function are the dimension of the Wiener
process, the step size and the desired precision:
```julia-repl
julia> alg = optimal_algorithm(m,h,err); # output: Fourier()
```

It is also possible to choose the algorithm directly using the
keyword argument `alg`. The value can be one of
`Fourier()`, `Milstein()`, `Wiktorsson()` and `MronRoe()`:
```julia-repl
julia> II = iterated_integrals(W,h; alg=Milstein())
```

As the norm for the considered error, e.g., the $\max,L^2$- and $\mathrm{F},L^2$-norm
can be selected using a keyword argument. The corresponding
values are `MaxL2()` and `FrobeniusL2()`:
```julia-repl
julia> II = iterated_integrals(W,h,err; error_norm=FrobeniusL2())
```

If iterated stochastic integrals for some $Q$-Wiener process need to
be simulated, like for the numerical simulation of solutions to SPDEs,
then the increment of the $Q$-Wiener process together with the
square roots of the eigenvalues of the associated covariance
operator have to be provided:
```julia-repl
julia> q = [1/k^2 for k=1:m]; # eigenvalues of cov. operator
julia> QW = sqrt(h) * sqrt.(q) .* randn(m); # Q-Wiener increment
julia> IIQ = iterated_integrals(QW,sqrt.(q),h,err)
```
In this case, the function `iterated_integrals` utilizes a
scaling of the iterated stochastic integrals and also adjusts the error
estimates appropriately such that the error bound holds w.r.t.\ the
iterated stochastic integrals $\mathcal{I}^{Q}(h)$ based on the
$Q$-Wiener process. Here the error norm defaults to the $\mathrm{F},L^2$-error.

Note that all discussed keyword arguments are optional and can be
combined as needed. 

Additional information can be found, e.g., using the Julia help mode:
```julia-repl
julia> ?iterated_integrals
julia> ?optimal_algorithm
```