using SRK
using Test
import Random
import Random: randn
import LinearAlgebra: diag

Random.seed!(638278)

@testset "Iterated Integrals" begin
    @testset "Expected Output m=5" begin
        W = randn(5)
        Ints = SRK.simdoubleintegrals_n(W,2) # 2 approximation terms
        expected_output = [0.418259  0.0677285 0.809848 -0.220337  -0.0419943;
                          -0.603668 -0.4218   -0.166696 -0.233224   0.0483685;
                           0.200639 -0.128188 -0.222006  0.0609429 -0.0259348;
                           0.546643  0.138     0.118597 -0.471011   0.0385036;
                          -0.583664  0.134213 -0.318314 -0.149669  -0.393427]
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
