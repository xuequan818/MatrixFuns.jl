const RealOrComplex = Union{Real,Complex}

include("split.jl")
include("reorder_schur.jl")
include("atomic_block_fun.jl")
include("parlett_recurrence.jl")

"""
    mat_fun(f, A; scale, sep, max_deg, color, ε, checknative)

Compute the matrix-variable function `f(A)`. Return a matrix.

`f` a general scalar function, and called as `f(x)`.

`A` an arbitrary square matrix.

`scale` the scaling of the Talyor series error, is used to control the spread of each splitting cluster.
By default, `scale = 1`.

`sep` the initial separation distance, to split the eigenvalues of `A` into clusters that satisfy
(i) min{|λ - μ|: λ ∈ cl_p, μ ∈ cl_q, p ≠ q} > sep.		
(ii) for cl_p with |cl_p| > 1, ∀ λ ∈ cl_p, ∃ μ ∈ cl_p and μ ≠ λ, s.t. |λ - μ| ≤ sep. 
Here cl_p is the cluster with index p.
By default, `sep = 0.1*scale`.

`max_deg` the maximum Taylor series order for the diagonal blocks computation.

`color(x::Number)::Integer` an index mapping function. 
This is to ensure all eigenvalues within a cluster have the same color index. 
By default, the `color` function maps `x` to `1`. 
For discontinuous functions, users can assign a distinct color to each continuous interval. 
E.g., for the `sign` function, users can customize the color mapping using 
`color = x -> Int(real(sign(x)))`

`ε` the target accuracy, set to machine accuracy by default.

`checknative` check `f` is a Julia native function.
"""
function mat_fun(f::Function, A::AbstractMatrix{TT};               
                 scale=1.0, sep=0.1scale, max_deg=250, 
                 color::Function=(x->1), 
                 ε=eps(real(float((TT)))), 
                 checknative=native(f)) where {TT<:RealOrComplex}      
    n = checksquare(A)

    if checknative
        fNM = f(NATIVE_TEST_MAT)
        if typeof(fNM) <: Number
            return Matrix(fNM*I,n,n)
        end
    end

    if isone(n)
        return elem_fun(f, A)
    end
    
    if ishermitian(A)
        if checknative
            return f(A)
        else
            Λ, Ψ = eigen(A)
            return diag_mat_fun(f, Λ, Ψ)
        end
    end
    
    # Schur decomposition of A
    sflag = 0
    S = schur(A)
    if !istriu(S.T)
        S = Schur{Complex}(schur(A))
        sflag = 1
    end
    T, Z, Λ = S

    # A is diagonalizable.
    if isdiag(T) || norm(T - Diagonal(Λ)) < ε/sep
        return diag_mat_fun(f, Λ, Z)
    end

    # split the eigenvalues into clusters
    split_map = get_splittings(Λ; sep, max_deg, scale, color, ε, checknative)

    # compute f(T)
    if maximum(split_map) == n
        F = parlett_recurrence(f, T)
    elseif allequal(split_map)
        F = atomic_block_fun(f, T; max_deg, checknative)
    else
        S, block = reorder_schur(S, split_map)
        F = block_parlett_recurrence(f, S.T, block; max_deg, tol_taylor=100ε, checknative)
    end

    FA = S.Z * F * S.Z'
    if isone(sflag)
        if TT <: Real && Base.return_types(f, Tuple{TT})[1] <: Real
            (norm(imag(FA)) < 100ε) && (FA = real(FA))
        end
    end

    return FA
end

diag_mat_fun(f::Function, Λ, P) = P * Diagonal(elem_fun(f, Λ)) * P'

# Check `f` is a Julia native function
const NATIVE_TEST_MAT = @MArray ones(1, 1)

function native(f::Function)
    try
        f(NATIVE_TEST_MAT)
        true
    catch
        false
    end
end
