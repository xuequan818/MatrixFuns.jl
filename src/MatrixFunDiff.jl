module MatrixFunDiff

using LinearAlgebra
using PermutationSymmetricTensors
using Combinatorics

export FiniteDual
export divided_difference
include("finitedual.jl")
include("divided_difference.jl")

export frechet_matrix_fun
include("matrixfun_derivative.jl")

end # module MatrixFunDiff
