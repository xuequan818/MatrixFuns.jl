function reorder_schur(S::Schur, split_map)
    rosp, reord, block = get_swappings(split_map)

	select = zeros(Bool, length(split_map))
	if !isnothing(reord) 
		for iro in reord
			@. select = 0
			@views irosp = rosp[iro]
			@. select[irosp] = 1
			ordschur!(S, select)
		end
	end

	return S, block
end

# Find the swap strategy that converts an unordered sequence             
# into a sequence that aligns the same indexes together.
# E.g., (1,4,2,1,2,3,3,2) -> (1,1,4,2,2,3,3,2) -> (2,2,2,1,1,4,3,3)			
function get_swappings(split_map::Vector{Int})
	N = length(split_map)

	# Filter out clusters with only one point and 
	# clusters that are already aligned together to reduce swaps.
    sp1 = sortperm(split_map)
	cl1(x) = searchsorted(split_map[sp1], x)
    ord = cl1.(unique(split_map))
	block_pos = [0]
	count1 = 0
	@views for io in ord
		isp = sp1[io]
		il = length(isp)
		if il == 1 || (maximum(isp)-minimum(isp)) == (il-1)
			count1 += il
			@. split_map[isp] = N+1
            push!(block_pos, count1)
		end
	end

    if count1 == N
        sp2, reord = nothing, nothing
	else
		# Find the swap strategy for the remaining clusters,
		# that converts an unordered sequence to a decreasing order sequence. 
		count2 = N - count1
		@. block_pos += count2
		sp2 = sortperm(split_map)
		sort_map2 = split_map[sp2]
		cl2(x) = searchsorted(sort_map2, x)
		reord = cl2.(unique(sort_map2[1:findlast(x -> x ≤ N, sort_map2)]))
		@views for iro in reord
			count2 -= length(iro)
			irosp = sp2[iro]
			for (i, il) in enumerate(irosp)
				irosp[i] += sum(x->x<split_map[il], split_map[il+1:end]; init=0)
			end
			pushfirst!(block_pos, count2)
		end
		@assert iszero(count2)
	end

    block = [block_pos[i]+1:block_pos[i+1] for i = 1:length(block_pos)-1]

    return sp2, reord, block
end
