module Swap

using Test
import MatrixFuns: get_swappings, split_cluster

const  N = [20, 30]

@testset "Swap strategy" begin
	for n in N
		asgmt = rand(1:10, n)
		reasgmt = copy(asgmt)
		resp, rerng, block = get_swappings(asgmt)

		if !isnothing(rerng)
			@views for irr in rerng
				irrsp = resp[irr]
				isp = vcat(irrsp, setdiff(1:n, irrsp))
				reasgmt = reasgmt[isp]
			end
		end

		@views for ib in block
			@test allequal(reasgmt[ib])
		end

		clsp, clrng = split_cluster(reasgmt)
		@test length.(clrng) == length.(block)
	end
end

end
