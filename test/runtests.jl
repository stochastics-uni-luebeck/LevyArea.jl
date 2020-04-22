using IteratedIntegrals
using Test
import Random
import Random: randn
import LinearAlgebra: diag

Random.seed!(638278)

@testset "Iterated Integrals" begin
#    @testset "Expected Output m=5" begin
#        W = randn(5)
#        Ints = IteratedIntegrals.simdoubleintegrals_n(W,2) # 2 approximation terms
#        expected_output = [0.418259  0.229801    1.42757  -0.374705  -0.304829
#                          -0.765741 -0.4218     -0.436137 -0.216972   0.180284
#                          -0.417082  0.141253   -0.222006 -0.133748   0.0398915
#                          0.701012   0.121748    0.313288 -0.471011  -0.0607856
#                          -0.320829  0.00229782 -0.38414  -0.0503793 -0.393427]
#        @test Ints ≈ expected_output atol=1e-5
#    end
    @testset "m=1, h=$h" for h in rand(5)
        W = √h * randn()
        @time Ints = simiterintegrals(W, h, h)
        @test Ints == 0.5W^2 - 0.5h
    end
    h = 0.01
    println("h = $h")
    @testset "Alg $alg: Diagonal m=$i" for alg in [Fourier(), Wiktorsson()], i in [1;50:50:500]
        W = √h * randn(i)
        print("Alg $alg: Diagonal m=$i ")
        @time Ints = simiterintegrals(W, h, h^(3/2), alg=alg) # stepsize h
        @test diag(Ints) ≈ 0.5*W.^2 .- 0.5h
    end
end
