module DividedDifference

using Test
using MatrixFuns
using DiffTests

# compute the divided difference by definition
function div_diff_def(f::Function, x::Vector{<:Real})
    n = length(x)
    isone(n) && return f(x[1])

    d0 = [(f(x[i]) - f(x[i+1])) / (x[i] - x[i+1]) for i = 1:n-1]
    for i = 2:n-1
        d1 = [(d0[j] - d0[j+1]) / (x[j] - x[j+i]) for j = 1:n-i]
        d0 = d1
    end
    @assert length(d0) == 1

    d0[1]
end

const N, X = 4, [1, 2, 3, 4]

@testset "$f" for f in DiffTests.NUMBER_TO_NUMBER_FUNCS
    for i = 1:N
		x = X[1:i]
		dd_def = div_diff_def(f, x)
		dd_table = div_diff(f, x)
        @test isapprox(dd_def, dd_table)
    end
end

@testset "1/(1+exp(1000*x))" begin
    f(x) = 1 / (1 + exp(1000 * (x-2.5)))
    color(x) = (x-2.5) < 0 ? 1 : 2
    for i = 1:N
        x = X[1:i]
        dd_def = div_diff_def(f, x)
        dd_table = div_diff(f, x; scale=0.01)
        @test isapprox(dd_def, dd_table)
    end
end

@testset "sign function" begin
    color(x) = Int(real(sign(x-2.5)))
    f(x) = sign(x-2.5)
    for i = 1:N
        x = X[1:i]
        dd_def = div_diff_def(f, x)
        dd_table = div_diff(f, x; sep=Inf, color)
        @test isapprox(dd_def, dd_table)
    end
end

@testset "heaviside step function" begin
    heaviside(x) = (x-2.5) < 0 ? 1.0 : 0.0
    color(x) = (x-2.5) < 0 ? 1 : 2
    for i = 1:N
        x = X[1:i]
        dd_def = div_diff_def(heaviside, x)
        dd_table = div_diff(heaviside, x; sep=Inf, color)
        @test isapprox(dd_def, dd_table)
    end
end

@testset "Arbitrary piecewise function" begin
    f(x) = x < 1.5 ? exp(1 / (x^2 + 1)) : (x < 2.5 ? sin(x) - 1 : 1)
    color(x) = x < 1.5 ? 1 : (x < 2.5 ? 2 : 3)
    for i = 1:N
        x = X[1:i]
        dd_def = div_diff_def(f, x)
        dd_table = div_diff(f, x; color)
        @test isapprox(dd_def, dd_table)
    end
end

@testset "complex output" begin
    f(x) = sin((1+im)*x)
    for i = 1:N
        i = 2
        x = X[1:i]
        dd_def = div_diff_def(f, x)
        dd_table = div_diff(f, x)
        @test isapprox(dd_def, dd_table)
    end
end

end
