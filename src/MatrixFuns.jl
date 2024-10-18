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

export mat_fun
include("matrix_function/matrix_function.jl")

end
