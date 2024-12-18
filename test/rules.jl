module MatFunChainRules

using Test
using MatrixFuns
using LinearAlgebra
using ChainRulesTestUtils

const N = 10

@testset "frule, T=$T" for T in (Float64, ComplexF64)
    A = Matrix(Hermitian(randn(T, N, N)))
    test_frule(mat_fun, x -> 1/(1+exp(x)) , A)
end

@testset "frule, complex out" begin
	A = Matrix(Hermitian(randn(ComplexF64, N, N)))
    test_frule(mat_fun, x -> sin((1+im)x), A)
end

@testset "rrule, T=$T" for T in (Float64, ComplexF64)
    A = Matrix(Hermitian(randn(T, N, N)))
    test_rrule(mat_fun, x -> 1 / (1 + exp(x)), A; check_inferred=false)
end

@testset "rrule, complex out" begin
    A = Matrix(Hermitian(randn(ComplexF64, N, N)))
    test_rrule(mat_fun, x -> sin((1 + im)x), A; check_inferred=false)
end

end