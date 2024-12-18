using MatrixFuns
using Test

@testset "MatrixFuns.jl" begin
    t0 = time()
    @testset "Split" begin
        println("##### Testing points splitting...")
        t = @elapsed include("split.jl")
        println("##### done (took $t seconds).")
    end
    @testset "Swap" begin
        println("##### Testing swap strategy ...")
        t = @elapsed include("swap.jl")
        println("##### done (took $t seconds).")
    end
    @testset "AtomicBlock" begin
        println("##### Testing atomic block computation functionality...")
        t = @elapsed include("atomic_block.jl")
        println("##### done (took $t seconds).")
    end
    @testset "MatrixFunction" begin
        println("##### Testing matrix-variable function computation functionality...")
        t = @elapsed include("matrix_function.jl")
        println("##### done (took $t seconds).")
    end
    @testset "DividedDifference" begin
        println("##### Testing divided difference functionality...")
        t = @elapsed include("divided_difference.jl")
        println("##### done (took $t seconds).")
    end
    @testset "DivDiffConvergence" begin
        println("##### Testing convergence of divided difference functionality...")
        t = @elapsed include("div_diff_convergence.jl")
        println("##### done (took $t seconds).")
    end
    @testset "Fréchet Derivative" begin
        println("##### Testing Fréchet derivative functionality...")
        t = @elapsed include("frechet_derivative.jl")
        println("##### done (took $t seconds).")
    end
    @testset "Chain Rules" begin
        println("##### Testing Chain Rules functionality...")
        t = @elapsed include("rules.jl")
        println("##### done (took $t seconds).")
    end
    println("##### Running all MatrixFuns tests took $(time() - t0) seconds.")
end
