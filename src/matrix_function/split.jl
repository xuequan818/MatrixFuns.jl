"""
    get_splittings(f, pts; sep, color, scale, max_deg, ε) -> Vector{Int}

Return the splitting cluster assignments vector z (z_i is the index of the splitting cluster for the i-th data points (eigenvalues)). 
"""
function get_splittings(f::Function, pts::Vector{T}; 
						sep=0.1, color::Function=(x->1),
                        checknative=native(f),
                        kwargs...) where {T<:Number}	
    # Split the points by `color`.
    asgmt = color.(pts)
    @assert eltype(asgmt) <: Integer
    clsp, clrng = split_cluster(asgmt)

    # Each point has a different color
    if isnothing(clrng)
        @. asgmt = 1:length(pts)
        return asgmt
    end

    if isinf(sep)
        @views for (i, ir) in enumerate(clrng)
            fill!(asgmt[clsp[ir]], i)
        end
        return asgmt
    end

    # Split each color cluster by separation 
    min_ind = 1
    @views for ir in clrng
        cl_i = clsp[ir]
        pts_i = pts[cl_i]
        asg_i = _get_splittings(pts_i, sep, min_ind, checknative; kwargs...)
        asgmt[cl_i] = asg_i
        min_ind = maximum(asg_i) + 1
    end

    return asgmt
end

# Repeatedly split the clusters by separation (halving the value at each iteration) 
# until the spread of subcluster meets the criteria 
function _get_splittings(pts::AbstractVector{T},
						 δ::Real, min_ind::Int,
                         checknative::Bool;
                         max_deg::Int=100,
						 scale::Real=1.0, 
						 ε=eps(real(T))) where {T<:Number}
    if length(pts) == 1
        return [min_ind]
    end

    asgmt, clsp, clrng = split_by_sep(pts, δ)
    @. asgmt += (min_ind - 1)

    if !isnothing(clrng)
        @views for ir in clrng
            cl_i = clsp[ir]
            pts_i = pts[cl_i]
            spread_i = get_spread(pts_i)
            can_split = checkspread(spread_i, scale, Val(checknative); max_deg, tol=ε/δ)
            asg_i = can_split ? _get_splittings(pts_i, δ/2, min_ind, checknative;
                				                max_deg, scale, ε) :
                    			asgmt[cl_i] .+ (min_ind - minimum(asgmt[cl_i]))
            asgmt[cl_i] = asg_i
            min_ind = maximum(asg_i) + 1
        end
    end

    return asgmt
end

# Split the same mapping index into one cluster
function split_cluster(asgmt::AbstractVector{T}; rngsort=true) where {T<:Integer}
    N = length(asgmt)

    if (minimum(asgmt) == 1 && maximum(asgmt) == N) || allunique(asgmt)
        return nothing, nothing
    end
	
    sp = sortperm(asgmt)
    asgmt_sort = asgmt[sp]
    searchclust(x) = searchsorted(asgmt_sort, x)
    rng = rngsort ? searchclust.(unique(asgmt_sort)) :
                    searchclust.(unique(asgmt))

    return sp, rng
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

    pos = vcat(0, findall(x -> x > δ, dist), N) # find all split positions
    asgmt = similar(pts, Int)
    nr = length(pos) - 1
    rng = Vector{UnitRange{Int64}}(undef, nr)
    l = 1
    for i = 1:nr
        ir = pos[i]+1:pos[i+1]
        rng[i] = ir
        fill!(asgmt[ir], l)
        l += 1
    end

    return asgmt[sortperm(sp)], sp, rng
end

function get_dist_vec(pts::AbstractVector{<:Real})
    # Sort the eigenvalues to quickly calculate distance
    sp = sortperm(pts; rev=true) # real pts
    pts_sp = pts[sp]
    @views dist = pts_sp[1:end-1] - pts_sp[2:end]

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
    asgmt = cutree(R, h=δ)
    sp, rng = split_cluster(asgmt)

    return asgmt, sp, rng
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

function checkspread(spread, scale, ::Val{false}; max_deg, tol)
    (spread / scale)^(max_deg + 1) > tol
end
checkspread(spread, scale, ::Val{true}; kwargs...) = spread > scale
