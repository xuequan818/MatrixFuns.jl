"""
Construct the splitting cluster assignments vector z 
(z_i is the index of the splitting cluster for the i-th data points (eigenvalues)). 
"""
function get_splittings(pts::Vector{T}; sep=0.1,
						color::Function=(x->1),
                        checknative=false,
                        kwargs...) where {T<:RealOrComplex}
    N = length(pts) 
    if N == 1
        return [1]
    end

    # Split the points by `color`.
    asgmt = color.(pts)
    @assert eltype(asgmt) <: Integer
    clsp, clrng = split_cluster(asgmt)

    # Each point has a different color
    if isnothing(clrng)
        @. asgmt = 1:N
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
                         max_deg::Int=250,
						 scale::Real=1.0, 
						 ε=eps(real(T))) where {T<:RealOrComplex}
    if length(pts) == 1
        return [min_ind]
    end

    asgmt, clsp, clrng = split_by_sep(pts, δ) # asgmt[clsp[clrng[i]]] .== i
    @. asgmt += (min_ind - 1)
    
    if isinf(scale)
        return asgmt
    end

    if !isnothing(clrng)
        @views for ir in clrng
            cl_i = clsp[ir]
            pts_i = pts[cl_i]
            can_split = checkspread(pts_i, scale, Val(checknative), max_deg, ε / δ)
            asg_i = can_split ? _get_splittings(pts_i, δ/2, min_ind, checknative;
                				                max_deg, scale, ε) :
                    			asgmt[cl_i] .+ (min_ind - minimum(asgmt[cl_i]))
            asgmt[cl_i] = asg_i
            min_ind = maximum(asg_i) + 1
        end
    end

    return asgmt
end


"""
Split the same color mapping index into one cluster
"""
function split_cluster(asgmt::AbstractVector{T}) where {T<:Integer}
    N = length(asgmt)

    if (minimum(asgmt) == 1 && maximum(asgmt) == N) || allunique(asgmt)
        return nothing, nothing
    end

    sp = issorted(asgmt) ? collect(1:N) : sortperm(asgmt)
    @views asgmt_sort = asgmt[sp]
    searchclust(x) = searchsorted(asgmt_sort, x)
    rng = searchclust.(unique(asgmt))         

    return sp, rng
end

"""
    split_by_sep(pts::AbstractVector{<:Real}, δ::Real)
    split_by_sep(pts::AbstractVector{<:Complex}, δ::Real)

Split the data points (eigenvalues) by separation parameter δ.
"""
function split_by_sep(pts::AbstractVector{<:Real}, δ::Real)
    N = length(pts)
    @assert N > 1

    # Sort the eigenvalues to quickly calculate distance
    dist, sp, decrease = get_dist_vec(pts)

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
    @views for i = 1:nr
        ir = pos[i]+1:pos[i+1]
        rng[i] = ir
        fill!(asgmt[ir], l)
        l += 1
    end
    decrease || (asgmt = asgmt[sortperm(sp)])

    return asgmt, sp, rng
end

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

function get_dist_vec(pts::AbstractVector{<:Real})
    # Sort points in decreasing order to quickly calculate distances
    if issorted(pts; rev=true)
        sp = collect(1:length(pts))
        decrease = true
    else
        sp = sortperm(pts; rev=true) 
        decrease = false
        pts = pts[sp]
    end
    @views dist = pts[1:end-1] - pts[2:end]

    return dist, sp, decrease
end

function get_dist_mat(px::AbstractVector{T}, 
                      py::AbstractVector{T}) where {T<:Number}
    map(x -> norm(x[1] - x[2]), Iterators.product(px, py))
end
get_dist_mat(pts::AbstractVector{<:Number}) = get_dist_mat(pts, pts)

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

"""
Find the spread of the data points (eigenvalues).
"""
get_spread(pts::AbstractVector{<:Real}) = maximum(pts) - minimum(pts)
get_spread(pts::AbstractVector{<:Complex}) = maximum(get_dist_mat(pts))
get_spread(pts::T) where {T<:Number} = zero(real(T))

"""
Check that the spread error is smaller than the splitting error.
"""
function checkspread(pts::AbstractVector, scale, isnative, max_deg, tol)
    spread = get_spread(pts)
    checkspread(spread, scale, isnative, max_deg, tol)
end
function checkspread(spread::Real, scale, ::Val{false}, max_deg, tol)
    (spread / scale)^(max_deg + 1) > tol
end
checkspread(spread::Real, scale, ::Val{true}, max_deg, tol) = spread > scale
