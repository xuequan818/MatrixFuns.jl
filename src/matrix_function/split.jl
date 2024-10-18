"""
    get_splittings(f, pts; sep, color, scale, max_deg, ε) -> Vector{Int}

Return the splitting cluster assignments vector z (z_i is the index of the splitting cluster for the i-th data points (eigenvalues)). 
"""
function get_splittings(f::Function, pts::Vector{T}; 
						sep=0.1, color::Function=(x->1),
                        isnative=native(f),
                        kwargs...) where {T<:Number}	
    # Split the points by `color`.
    clsp, clrg = split_cluster(color.(pts))

	if isnothing(clrg)
        # Each point has a different color 
		split_map = collect(1:length(pts))
	else
        # Split each cluster by separation 
        min_map = 1
		split_map = similar(pts, Int)
        @views for ir in clrg
            cl_i = clsp[ir]
			pts_i = pts[cl_i]
            map_i = _get_splittings(pts_i, sep, min_map, isnative; kwargs...)
			split_map[cl_i] = map_i
            min_map = maximum(map_i) + 1
        end
	end

	split_map
end

# Repeatedly split the clusters by separation (halving the value at each iteration) 
# until the spread of subcluster meets the criteria 
function _get_splittings(pts::AbstractVector{T},
						 δ::Real, min_map::Int,
                         isnative::Bool;
                         max_deg::Int=100,
						 scale::Real=1.0, 
						 ε=eps(real(T))) where {T<:Number}
    if length(pts) == 1
        return [min_map]
    end

    split_map, clsp, clrg = split_by_sep(pts, δ)
    @. split_map += (min_map - 1)

    if !isnothing(clrg)
        @views for ir in clrg
            cl_i = clsp[ir]
            pts_i = pts[cl_i]
            spread_i = get_spread(pts_i)

            splitting = false
            if isnative
                (spread_i > scale) && (splitting = true)
            elseif (spread_i/scale)^(max_deg+1) > ε/δ
                splitting = true
            end

            map_i = splitting ? _get_splittings(pts_i, δ/2, min_map, isnative;
                				max_deg, scale, ε) :
                    			split_map[cl_i] .+ (min_map - minimum(split_map[cl_i]))
            split_map[cl_i] = map_i
            min_map = maximum(map_i) + 1
        end
    end

    return split_map
end

# Split the same mapping index into one cluster
function split_cluster(index::AbstractVector{T}) where {T<:Integer}
    N = length(index)

    if (minimum(index) == 1 && maximum(index) == N) || unique(index) == index
        return nothing, nothing
    end
	
    sp = sortperm(index)
    index_sort = index[sp]
    cl(x) = searchsorted(index_sort, x)

    return sp, cl.(unique(index_sort))
end

# Real points
function split_by_sep(pts::AbstractVector{<:Real}, δ::Real)
    N = length(pts)
    @assert N > 1

    # Sort the eigenvalues to quickly calculate distance
    dist, sp = get_dist_vec(pts)

    if minimum(dist) > δ # All points are separated.
        return collect(1:N), nothing, nothing
	end

    if maximum(dist) ≤ δ # All points are in one cluster.
        return ones(Int, N), sp, [1:N]
	end 

    split_pos = vcat(0, findall(x -> x > δ, dist), N) # find all split positions
    split_map = similar(pts, Int)
    nr = length(split_pos) - 1
    split_range = Vector{UnitRange{Int64}}(undef, nr)
    l = 1
    for i = 1:nr
        ir = split_pos[i]+1:split_pos[i+1]
        split_range[i] = ir
        @. split_map[ir] = l
        l += 1
    end

    return split_map[sortperm(sp)], sp, split_range
end

function get_dist_vec(pts::AbstractVector{<:Real})
    # Sort the eigenvalues to quickly calculate distance
    sp = sortperm(pts; rev=true) # real pts
    pts_sp = pts[sp]
    dist = pts_sp[1:end-1] - pts_sp[2:end]

    return dist, sp
end

# Complex points
function split_by_sep(pts::AbstractVector{<:Complex}, δ::Real)
    N = length(pts)
    @assert length(pts) > 1

    dist = get_dist_mat(pts)

    max_gap = maximum(dist)
    if max_gap ≤ δ # All points are in one cluster.
        return ones(Int, N), collect(1:N), [1:N]
    end

    retriu!(dist, 1; δ=max_gap)
    if minimum(dist) > δ # All points are separated.
        return collect(1:N), nothing, nothing
    end

    R = hclust(dist, linkage=:single, uplo=:U)
    split_map = cutree(R, h=δ)
    sp, split_range = split_cluster(split_map)

    return split_map, sp, split_range
end

function get_dist_mat(pts::AbstractVector{<:Number})
    map(x -> norm(x[1] - x[2]), Iterators.product(pts, pts))
end

# Keep the upper triangle of a matrix and 
# let the other parts be `δ`, overwriting `M` in the process.
function retriu!(M::AbstractMatrix{T}, k::Integer; δ=zero(T)) where {T}
    require_one_based_indexing(M)
    m, n = size(M)
    for j in 1:min(n, m + k)
        for i in max(1, j - k + 1):m
            M[i, j] = δ
        end
    end
    M
end

get_spread(pts::AbstractVector{<:Real}) = maximum(pts) - minimum(pts)
get_spread(pts::AbstractVector{<:Complex}) = maximum(get_dist_mat(pts))
get_spread(pts::T) where {T<:Number} = zero(real(T))
