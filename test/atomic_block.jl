module AtomicBlock

using Test
using LinearAlgebra
import MatrixFuns: atomic_block_fun

const N, NATIVE_FUNS = [20, 50], [exp, sin]

@testset "$f" for f in NATIVE_FUNS
	for n in N
		A = triu!(rand(n,n))
		ref = f(A)
		result = atomic_block_fun(f, A; checknative=false)
		@test isapprox(ref, result)
	end		
end

@testset "1/(1+exp(x))" begin
    f(x) = 1 / (1 + exp(x))
	for n in N
		A = triu!(rand(n,n))
	    ref = I / (I + exp(A))
		result = atomic_block_fun(f, A)
   	 	@test isapprox(ref, result)
	end
end

@testset "complex output" begin
    f(x) = sin((1 + im) * x)
    for n in N
        A = triu!(rand(n, n))
        ref = f(A)
        result = atomic_block_fun(f, A)
        @test isapprox(ref, result)
    end
end

end
