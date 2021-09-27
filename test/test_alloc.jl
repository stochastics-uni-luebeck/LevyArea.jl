using IteratedIntegrals
import IteratedIntegrals: levyarea
import Random: default_rng, randn!
using SpecialFunctions: trigamma

function levyarea2(W::AbstractVector{T}, n::Integer, alg::Fourier) where {T<:AbstractFloat}
    rng = default_rng()
    m = length(W)
    X = randn(rng, T, n, m) # allocates m*n Floats
    Y = randn(rng, T, m, n) # allocates m*n Floats
    Y .= (Y .- √(T(2)).*W) ./ (1:n)'
    A = Y * X # allocates m*m Floats
    # Antisymmetrize
    # G = inv(2*T(π)).*(A .- A') # allocates m*m Floats
    for i = 1:m
        @inbounds A[i,i] = zero(T)
        for j = (i+1):m
            @inbounds A[i,j] = (A[i,j] - A[j,i])/(2*T(π))
            @inbounds A[j,i] = -A[i,j]
        end
    end
    # @assert A ≈ G
    return A
end
function levyarea2(W::AbstractVector{T}, n::Integer, alg::Milstein) where {T<:AbstractFloat}
    rng = default_rng()
    m = length(W)
    X = randn(rng, T, n, m) # allocates m*n Floats
    Y = randn(rng, T, m, n) # allocates m*n Floats
    Y .= (Y .- √(T(2)).*W) ./ (1:n)'
    A = Y * X # allocates m*m Floats
    # Add simple rest approximation (a₀)
    # M = randn(rng, T, m) # allocates m Floats
    M = randn!(rng, view(Y, :, 1))
    a = T(√(2*trigamma(n+1)))
    A .+= a .* W .* M'
    # Antisymmetrize
    # G = inv(2*T(π)).*(A .- A') # allocates m*m Floats
    for i = 1:m
        @inbounds A[i,i] = zero(T)
        for j = (i+1):m
            @inbounds A[i,j] = (A[i,j] - A[j,i])/(2*T(π))
            @inbounds A[j,i] = -A[i,j]
        end
    end
    # @show A ≈ G
    return A
end

T = Float16
W = randn(T, 100)
n = 1000

for alg ∈ [Fourier(), Milstein()]
    @show alg
    @show typeof(levyarea(W, n, alg)) == Matrix{T}
    # @code_warntype levyarea(1.0:3.0, n, Wiktorsson())
    # @show @ballocated levyarea(W, n, $alg)
    # @show @belapsed levyarea(W, n, $alg)
    b = @benchmark levyarea($W, $n, $alg)
    @show minimum(b)
    @show b.memory
end
## m=5, n=100
# 32.58, 33.07, 33.20, 32.78 μs
# 2800, 2896, 2896, 2656 bytes
## m=100, n=1000
# 48.986, 49.103, 49.100, 49.094 ms
# 440720, 441008, 441008, 420608
for alg ∈ [Fourier(), Milstein()]
    @show alg
    @show typeof(levyarea2(W, n, alg)) == Matrix{T}
    b = @benchmark levyarea2($W, $n, $alg)
    @show minimum(b)
    @show b.memory
end
## m=5, n=100
# 32.32, 32.77 μs
# 2656, 2752 bytes
## m=100, n=1000
# 49.010, 48.921 ms
# 420608, 420896
