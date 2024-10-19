module MatrixFuns

using StaticArrays
using LinearAlgebra
using LinearAlgebra: checksquare
using Base: require_one_based_indexing
using IterTools
using Clustering
using TaylorSeries
using Arblib
using SpecialFunctions
using Combinatorics
using PermutationSymmetricTensors

export mat_fun
include("matrix_function/matrix_function.jl")

export div_diff
export mat_fun_frechet
include("divided_difference.jl")
include("frechet_derivative.jl")

end
