module DivDiffConvergence

using Test
using MatrixFuns
using ForwardDiff
using DiffTests

# compute the n-th order derivative
function nth_derivative(f::Function, x::Number, n::Integer)
    iszero(n) ? f(x) : ForwardDiff.derivative(x -> nth_derivative(f, x, n-1), x)
end

const x, N = 0.1, 4

@testset "$f" for f in DiffTests.NUMBER_TO_NUMBER_FUNCS
    for i = 1:N
        d = nth_derivative(f, x, i - 1)
        dd = div_diff(f, fill(x, i))
        @test isapprox(d / factorial(i - 1), dd)
    end
end

@testset "1/(1+exp(1000*x))" begin
    f(x) = 1 / (1 + exp(1000 * x))
    color(x) = x < 0 ? 1 : 2
    for i = 1:N
        d = nth_derivative(f, x, i - 1)
        dd = div_diff(f, fill(x, i); color)
        @test isapprox(d / factorial(i - 1), dd)
    end
end

@testset "sign function" begin
    color(x) = Int(real(sign(x)))
    for i = 1:N
        d = nth_derivative(sign, x, i - 1)
        dd = div_diff(sign, fill(x, i); color)
        @test isapprox(d / factorial(i - 1), dd)
    end
end

@testset "heaviside step function" begin
    heaviside(x) = x < 0 ? 0 : 1
    color(x) = x < 0 ? 1 : 2
    for i = 1:N
        d = nth_derivative(heaviside, x, i - 1)
        dd = div_diff(heaviside, fill(x, i); color)
        @test isapprox(d / factorial(i - 1), dd)
    end
end

@testset "Arbitrary piecewise function" begin
    f(x) = x < 1.5 ? exp(1 / (x^2 + 1)) : (x < 2.5 ? sin(x) - 1 : 1)
    color(x) = x < 1.5 ? 1 : (x < 2.5 ? 2 : 3)
    for i = 1:N
        d = nth_derivative(f, x, i - 1)
        dd = div_diff(f, fill(x, i); color)
        @test isapprox(d / factorial(i - 1), dd)
    end
end

@testset "complex output" begin
    f(x) = sin((1 + im) * x)
    for i = 1:N
        d = nth_derivative(f, x, i - 1)
        dd = div_diff(f, fill(x, i))
        @test isapprox(d / factorial(i - 1), dd)
    end
end

end
