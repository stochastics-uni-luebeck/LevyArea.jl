using SRK
using Test
import Random
import Random: randn
import LinearAlgebra: diag

Random.seed!(638278)

@testset "Iterated Integrals" begin
    @testset "Expected Output m=5" begin
        W = randn(5)
        Ints = simdoubleintegrals(W,2)
        expected_output = [ 0.418259   -0.560299    0.14876    0.390997   -1.33136;
                            0.0243588  -0.4218     -0.208187  -0.334559    0.357775;
                            0.861727   -0.0866971  -0.222006  -0.0437886  -0.224937;
                           -0.0646903   0.239335    0.223329  -0.471011   -0.258644;
                            0.705706   -0.175193   -0.119312   0.147479   -0.393427]
        @test Ints ≈ expected_output atol=1e-5
    end
    @testset "Diagonal m=$i" for i in 5:5:50
        W = randn(i)
        Ints = simdoubleintegrals(W, 1.0, 1.0)
        @test diag(Ints) == 0.5*W.^2 .- 0.5
    end
end
