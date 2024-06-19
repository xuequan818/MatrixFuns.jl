module MatrixFunDiff

using ForwardDiff
using LinearAlgebra

export FiniteDual
export divided_difference
include("finitedual.jl")
include("divided_difference.jl")

include("matrix_derivative.jl")

end # module MatrixFunDiff
