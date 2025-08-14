using BenchmarkTools
using LinearAlgebra
using SpecialFunctions
using MatrixFuns

name(f) = last(split(string(f), '.'))

const SUITE = BenchmarkGroup()

inv_exp(x) = 1 / (1 + exp(100 * x))
const funcs = [exp, erf, inv_exp]
const scales = [1, 1, 1/100]
const sizes = [5, 50, 100]
const A = [rand(n, n) for n in sizes]
const X = [rand(n) for n in sizes]

const H = map([5, 20, 80]) do n; X = rand(ComplexF64, n, n); H = 0.5 * (X + X') end
const hs = [map(i -> i * h, [1, 2]) for h in H]

const matrix_function_group = addgroup!(SUITE, "Matrix Function")
const divided_difference_group = addgroup!(SUITE, "Divided Difference")
const frechet_derivative_group = addgroup!(SUITE, "Fréchet Derivative")


for (f, s) in zip(funcs, scales)
	fm = addgroup!(matrix_function_group, name(f))
	dd = addgroup!(divided_difference_group, name(f))
	fd = addgroup!(frechet_derivative_group, name(f))
	for (i, n) in enumerate(sizes)
		Ai = A[i]
        fm[n] = @benchmarkable mat_fun($f, $Ai; scale=$s)

		Xi = X[i]
        dd[n] = @benchmarkable div_diff($f, $Xi; scale=$s)
	end

	for (i, Hi) in enumerate(H)
        hsi = hs[i]
		n = size(Hi, 1)
        fd[n] = @benchmarkable mat_fun_frechet($f, $Hi, $hsi; scale=$s)
	end
end
