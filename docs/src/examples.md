# Examples

## A simple Milstein scheme

We consider the following SDE
```math
\begin{aligned}
\mathrm{d}X^1_t &= \mathrm{d}W^1_t \\
\mathrm{d}X^2_t &= \mathrm{d}W^2_t \\
\mathrm{d}X^3_t &= X^1\mathrm{d}W^2_t\\
X_0 &= 0
\end{aligned}
```
with $X_t\in\mathbb{R}^3$, $t\in[0,1]$ and a two-dimensional Wiener process $W_t$.
This is a simple example of a non-commutative SDE.
The exact solution is
```math
\begin{aligned}
X^1_t &= W^1_t \\
X^2_t &= W^2_t \\
X^3_t &= \int_{0}^{t}\int_{0}^s\mathrm{d}W^1_r\mathrm{d}W^2_s \,.
\end{aligned}
```
We can see that the third component computes the Lévy area defined by the Wiener process.

Using the Euler-Maruyama scheme we can approximate this solution as
```math
\begin{aligned}
Y^1_{n+1} &= Y^1_{n} + \Delta W^1_n \\
Y^2_{n+1} &= Y^2_{n} + \Delta W^2_n \\
Y^3_{n+1} &= Y^3_{n} + Y^1_{n} \cdot \Delta W^2_n
\end{aligned}
```
whereas the Milstein scheme takes the following form
```math
\begin{aligned}
Z^1_{n+1} &= Z^1_{n} + \Delta W^1_n \\
Z^2_{n+1} &= Z^2_{n} + \Delta W^2_n \\
Z^3_{n+1} &= Z^3_{n} + Z^1_{n} \cdot \Delta W^2_n + \int_{t_n}^{t_{n+1}}\int_{t_n}^s\mathrm{d}W^1_r\mathrm{d}W^2_s \,.
\end{aligned}
```

The following example shows how to implement these two schemes.
We determine the optimal algorithm first, so it does not have to be recomputed every iteration.
```@example
using LevyArea
using Plots

T = 1
N = 100
Δt = T/N
Y₀ = Z₀ = [0.0,0.0,0.0]

Y = zeros(3,N+1)
Z = zeros(3,N+1)
Y[:,1] .= Y₀
Z[:,1] .= Z₀

# which algorithm is optimal for a two-dimensional
# Wiener process and the chosen stepsize
alg = optimal_algorithm(2, Δt)

for n = 1:N
    # simulate increment and iterated integral
    ΔW = sqrt(Δt)*randn(2)
    II = iterated_integrals(ΔW, Δt; alg=alg)

    # Euler-Maruyama scheme
    Y[1,n+1] = Y[1,n] + ΔW[1]
    Y[2,n+1] = Y[2,n] + ΔW[2]
    Y[3,n+1] = Y[3,n] + Y[1,n]*ΔW[2]

    # Milstein scheme
    Z[1,n+1] = Z[1,n] + ΔW[1]
    Z[2,n+1] = Z[2,n] + ΔW[2]
    Z[3,n+1] = Z[3,n] + Z[1,n]*ΔW[2] + II[1,2]
end

p_WP = plot(0:Δt:T, Y[1:2,:]', label=["\$X^1\$" "\$X^2\$"], title="Wiener process", legend=:topleft)
p_LA = plot(0:Δt:T, Y[3,:], label="Euler-Maruyama", title="Lévy area", legend=:topleft)
plot!(p_LA, 0:Δt:T, Z[3,:], label="Milstein")

plot(p_WP, p_LA, layout=@layout [a; b])
```
The upper plot shows the first two components of the solution, which are the same for both schemes.
The lower plot then shows the different approximations of the Lévy area (third component of the solution).