function simdoubleintegrals_old(W::AbstractVector{<:AbstractFloat}, n::Integer)
    # Preparation
    m::Integer = length(W)
    M::Integer = m*(m-1) ÷ 2
    Pₘ = sparse(1:m^2, cld.(1:m^2, m).+m.*mod.(0:(m^2-1), m), 1)
    Kₘ = sparse(1:M, [i+(j-1)*m for i=1:m, j=1:m if j<i], 1, M, m^2)
    KIP = Kₘ*(I-Pₘ)

    # 1. Simulate Brownian motion and determine n
    # These are given as input arguments.
    # 2.a Simulate Xₖ and Yₖ
    X = randn(m,n)
    YW = randn(m,n) .+ √(2).*W
    # 2.b Approximate stochastic area integral
    Aₙ = X*Diagonal(1 ./ (1:n))*YW'
    Aₙ = Aₙ'-Aₙ
    # 3. Simulate Gₙ and add the tail-sum approximation
    temp = KIP * kron(sparse(1.0I,m,m),W)
    Σ = 2.0 * temp * temp' + 2.0I
    sqrtΣ = (Σ + (2*√(1+W'*W))*I)
    lmul!(1 / (√2 + √(2+2*W'*W)), sqrtΣ)
    Aₙ .= Aₙ .+ √trigamma(n+1) .* reshape( KIP'*sqrtΣ*randn(M) ,m,m)
    # 4. Calculate the iterated integrals
    Iᵢⱼ = (W*W'-I) ./2 .+ 1/2π .* Aₙ
    # display(Aₙ)
    return Iᵢⱼ
end

function simdoubleintegrals_old2(W::AbstractVector{<:AbstractFloat}, n::Integer;
    rng = GLOBAL_RNG)
    m::Integer = length(W)
    # 1. Simulate Brownian motion and determine n
    #    These are given as input arguments.
    # 2. Simulate Xₖ and Yₖ and approximate stochastic area integral
    A = (randn(rng,m,n) .+ √(2).*W) ./ (1:n)' * randn(rng,n,m)
    # 3.a Simulate Gₙ
    G = zeros(m,m)
    G[tril!(trues(m,m),-1)] = √(2*trigamma(n+1)) * randn(rng,m*(m-1)÷2)
    # 3.b and add the tail-sum approximation
    A .+= 1 ./ (1+√(1+W'*W)) .* W*W'*(G-G') .+ G
    A = A - A'
    # 4. Calculate the iterated integrals
    (W*W'-I)./2 .+ A./2π
end
