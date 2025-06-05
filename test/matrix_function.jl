module MatrixFunction

using Test
using LinearAlgebra
using MatrixFuns

const N, SEP, NATIVE_FUNS = [20, 50], [0.1, 0.05], [exp, sin]

@testset "$f" for f in NATIVE_FUNS
    for (n, sep) in zip(N, SEP)
        A = rand(n,n)
        ref = f(A)
        result = mat_fun(f, A; sep, checknative=false)
        @test isapprox(ref, result)
    end
end

@testset "1/(1+exp(1000*x))" begin
    f(x) = 1 / (1 + exp(1000 * x))
    color(x) = x < 0 ? 1 : 2
    for n in N
        D = rand(n) .- 0.5
        (0 in D) && (D[findall(iszero, D)] .= 0.1)
        P = rand(n, n)
        A = P * Diagonal(D) * inv(P)
        ref = P * Diagonal(f.(D)) * inv(P)
        result = mat_fun(f, A; scale=0.001, checknative=false, color)
        @test isapprox(ref, result)
    end
end

@testset "sign function" begin
    color(x) = Int(real(sign(x)))
	for n in N
		D = rand(n) .- 0.5
		if 0 in D 
			D[findall(iszero, D)] .= 0.5
		end
		P = rand(n,n)
		A = P * Diagonal(D) * inv(P)
		ref = P * Diagonal(sign.(D)) * inv(P)
		result = mat_fun(sign, A; color)
		@test isapprox(ref, result)
	end
end

@testset "heaviside step function" begin
    heaviside(x) = x < 0 ? 0 : 1
	color(x) = x < 0 ? 1 : 2
    for n in N
        D = rand(n) .- 0.5
        if 0 in D
            D[findall(iszero, D)] .= 0.5
        end
        P = rand(n, n)
        A = P * Diagonal(D) * inv(P)
        ref = P * Diagonal(heaviside.(D)) * inv(P)
        result = mat_fun(heaviside, A; color)      
        @test isapprox(ref, result)
    end
end

@testset "Arbitrary piecewise function" begin
    f(x) = x < 0 ? 1/(x-1) : (x < 0.5 ? sin(x) : 1)
    color(x) = x < 0 ? 1 : (x < 0.5 ? 2 : 3)
    for (n, sep) in zip(N, SEP)
        D = 2 .* (rand(n) .- 0.5)
        (0 in D) && (D[findall(iszero, D)] .= 0.1)
        (0.5 in D) && (D[findall(isequal(0.5), D)] .= 0.6)
        P = rand(n, n)
        A = P * Diagonal(D) * inv(P)
        ref = P * Diagonal(f.(D)) * inv(P)
        result = mat_fun(f, A; sep, color)
        @test isapprox(ref, result)
    end
end

@testset "complex output" begin
    f(x) = sin((1 + im) * x)
    for n in N
        A = rand(n, n)
        ref = f(A)
        result = mat_fun(f, A)
        @test isapprox(ref, result)
    end
end

end