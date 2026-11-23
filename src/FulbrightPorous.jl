module FulbrightPorous
using Reexport

@reexport using LyoPronto
import LyoPronto: ODEProblem
@reexport using DrWatson
@reexport using Unitful
@reexport using Plots
@reexport using StatsPlots: @df, groupedbar, groupedbar!
@reexport using LaTeXStrings
@reexport using CSV
@reexport using TypedTables
@reexport using PrettyTables
@reexport using Optim
@reexport using NonlinearSolve
@reexport using OptimizationOptimJL
@reexport using Accessors
@reexport using DataInterpolations
@reexport using Latexify
@reexport using SavitzkyGolay
@reexport using TransformVariables
using ADTypes: AutoForwardDiff
using DocStringExtensions
using ConcreteStructs
using LinearAlgebra: Diagonal
using NaNMath
using StaticArrays: SA, SArray

export OP, optalg
const OP = OptimizationProblem
using LineSearches
optalg = BFGS(linesearch=LineSearches.BackTracking())

export odealg_chunk4
const odealg_chunk4 = LyoPronto.get_odealg(4)


include("pikal_microcollapse.jl")
include("kinetic_models.jl")
include("fitting.jl")
include("table_renaming.jl")
export blankplot_hrC, blankplot_lRp, fp_plot_defaults 
include("recipes.jl")

end