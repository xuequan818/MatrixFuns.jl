function reorder_schur(S::Schur, asgmt)
    resp, rerng, block = get_swappings(asgmt)

	select = zeros(Bool, length(asgmt))
	if !isnothing(rerng) 
		@views for irr in rerng
            fill!(select, 0)
			irrsp = resp[irr]
			fill!(select[irrsp], 1)
			ordschur!(S, select)
		end
	end

	return S, block
end

# Find the swap strategy that converts an unordered sequence             
# into a sequence that aligns the same indexes together.
# E.g., (1,4,2,1,2,3,3,2) -> (1,1,4,2,2,3,3,2) -> (2,2,2,1,1,4,3,3)			
function get_swappings(asgmt::Vector{Int})
	N = length(asgmt)

	# Filter out clusters with only one point and 
	# clusters that are already aligned together to reduce swaps.
	clsp, clrng = split_cluster(asgmt; rngsort=false)
	pos = [0]
	count1 = 0
	@views for ir in clrng
        cl_i = clsp[ir]
        li = length(cl_i)
        if li == 1 || (maximum(cl_i) - minimum(cl_i)) == (li- 1)
			count1 += li
            fill!(asgmt[cl_i], N+1)
            push!(pos, count1)
		end
	end

    # Find the swap strategy for the remaining clusters,
    # that converts an unordered sequence to a decreasing order sequence. 
    if count1 == N
        resp, rerng = nothing, nothing
	else
		count2 = N - count1
		@. pos += count2
		resp = sortperm(asgmt)
        reasgmt = asgmt[resp]
        searchreord(x) = searchsorted(reasgmt, x)
        rerng = searchreord.(unique(reasgmt[1:findlast(x -> x ≤ N, reasgmt)]))
        @views for irr in rerng
			count2 -= length(irr)
			irrsp = resp[irr]
            for (j, jr) in enumerate(irrsp)
				irrsp[j] += sum(x->x<asgmt[jr], asgmt[jr+1:end]; init=0)
			end
			pushfirst!(pos, count2)
		end
		@assert iszero(count2)
	end

    block = [pos[i]+1:pos[i+1] for i = 1:length(pos)-1]

    return resp, rerng, block
end
