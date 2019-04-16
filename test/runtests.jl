using IteratedIntegrals
using Test
import Random
import Random: randn
import LinearAlgebra: diag

Random.seed!(638278)

@testset "Iterated Integrals" begin
    @testset "Expected Output m=5" begin
        W = randn(5)
        Ints = IteratedIntegrals.simdoubleintegrals_n(W,2) # 2 approximation terms
        expected_output = [0.418259  0.229801    1.42757  -0.374705  -0.304829
                          -0.765741 -0.4218     -0.436137 -0.216972   0.180284
                          -0.417082  0.141253   -0.222006 -0.133748   0.0398915
                          0.701012   0.121748    0.313288 -0.471011  -0.0607856
                          -0.320829  0.00229782 -0.38414  -0.0503793 -0.393427]
        @test Ints ≈ expected_output atol=1e-5
    end
    @testset "m=1, h=$h" for h in rand(5)
        W = √h * randn()
        @time Ints = simdoubleintegrals(W,h)
        @test Ints == 0.5W^2 - 0.5h
    end
    @testset "Diagonal m=$i, h=$h" for i in [1;50:50:500], h in rand(2)
        W = √h * randn(i)
        @time Ints = simdoubleintegrals(W, h) # stepsize h
        @test diag(Ints) ≈ 0.5*W.^2 .- 0.5h
    end
end
