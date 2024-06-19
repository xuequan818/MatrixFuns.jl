using Test
using MatrixFunDiff

function dd_naive(f::Function, x...)
    N = length(x)

    d1 = [(f(x[i]) - f(x[i+1])) / (x[i] - x[i+1]) for i = 1:N-1]
    dd = Vector[d1]
    for i = 2:N-1
        Ni = N - i
        @views di = dd[i-1]
        push!(dd, [(di[j] - di[j+1]) / (x[j] - x[j+i]) for j = 1:Ni])
    end

    dd[end][1]
end

f(x) = x <= 0 ? 1/(1+exp(20*sin(x^2))) : 0
xx = collect(-1:3)
for i = 1:length(xx)-1
	for j = 1 : length(xx)-i
        @views xij = xx[j:j+i]
		@test abs(dd_naive(f,xij...)-divided_difference(f,xij...)) < 1e-10
	end
end