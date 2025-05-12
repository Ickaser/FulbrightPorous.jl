module FulbrightPores
using Reexport

@reexport using LyoPronto
import LyoPronto: ODEProblem
@reexport using DrWatson
@reexport using Unitful
@reexport using Plots
@reexport using StatsPlots: @df
@reexport using LaTeXStrings
@reexport using CSV
@reexport using TypedTables
@reexport using PrettyTables
@reexport using Optim
# @reexport using LeastSquaresOptim
@reexport using NonlinearSolve
@reexport using Optimization
@reexport using OptimizationOptimJL
@reexport using Accessors
@reexport using DataInterpolations
@reexport using UnitfulLatexify, Latexify
@reexport using SavitzkyGolay
using ADTypes: AutoForwardDiff
using DocStringExtensions
using LinearAlgebra: Diagonal


export odealg_chunk4
const odealg_chunk4 = Rodas4(autodiff=AutoForwardDiff(chunksize=4))

include("pikal_microcollapse.jl")
include("kinetic_models.jl")
include("fitting.jl")
include("table_renaming.jl")
include("recipes.jl")

end