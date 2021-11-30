using IteratedIntegrals
using Test, BenchmarkTools
import Random
import LinearAlgebra: diag

BenchmarkTools.DEFAULT_PARAMETERS.seconds = 1
Random.seed!(638278)

@testset "Iterated Integrals" begin
    @testset "m=1, h=$h" for h in rand(5)
        W = √h * randn()
        Ints = iterated_integrals(W, h, h^(3/2))
        @test Ints == 0.5W^2 - 0.5h
    end

    #trigger compilation
    for alg ∈ IteratedIntegrals.ITER_INT_ALGS
        iterated_integrals(randn(2),0.1,0.1, alg=alg)
    end
    
    h = 0.0001
    ϵ = h^(3/2)
    println("h = $h, ϵ = $ϵ")
    @testset "Alg $alg: m=$m" for alg ∈ IteratedIntegrals.ITER_INT_ALGS, m in [2;10;50;100;500]
        W = √h * randn(m)
        print("Alg $alg: m=$m n=$(terms_needed(m,h,ϵ,alg,MaxL2())) cost=$(IteratedIntegrals.effective_cost(m,h,ϵ,alg,MaxL2()))")
        @btime iterated_integrals($W, $h, $ϵ, alg=$alg) # stepsize h
        Ints = iterated_integrals(W, h, ϵ, alg=alg) # stepsize h
        @test diag(Ints) ≈ 0.5*W.^2 .- 0.5h
    end
end
