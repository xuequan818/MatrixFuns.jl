module Split

using Test
using LinearAlgebra
import MatrixFuns: split_by_sep, get_splittings, get_dist_mat, get_spread

const SEP, SCALE, TT = [0.05,0.1], [0.1,1], [Float64, ComplexF64]

@testset "Split by separation" begin
    for sep in SEP, T in TT
        pts = rand(T, 20)
        asgmt, _, _ = split_by_sep(pts, sep)
        bks = [findall(isequal(i), asgmt) for i = 1:maximum(asgmt)]
        @views for (i, ibk) in enumerate(bks)
            ipts = pts[ibk]
            if length(ipts) > 1
				flag = 0
				for (k, pk) in enumerate(ipts)
					for (l, pl) in enumerate(ipts)
						(k != l && norm(pk-pl) ≤ sep) && (flag = 1)
					end
				end
				@test isone(flag)
            end			
            for (j, jbk) in enumerate(bks)
                if j != i
                    jpts = pts[jbk]
                    dmatij = get_dist_mat(ipts, jpts)
                    @test minimum(dmatij) > sep
                end
            end
        end
    end
end

@testset "Split by separation and spread" begin
	for sep in SEP, scale in SCALE, T in TT
		pts = rand(T, 20)
		asgmt = get_splittings(pts; sep, scale, checknative=true)
		bks = [findall(isequal(i), asgmt) for i = 1:maximum(asgmt)]
		@views for (i, ibk) in enumerate(bks)
			ipts = pts[ibk]
            if length(ipts) > 1
                flag = 0
                for (k, pk) in enumerate(ipts)
                    for (l, pl) in enumerate(ipts)
                        (k != l && norm(pk-pl) ≤ sep) && (flag = 1)
                    end
                end
                @test isone(flag)
            end
			@test get_spread(ipts) ≤ scale
		end
	end
end

end
