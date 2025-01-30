@quickactivate :FulbrightPores

const OP = OptimizationProblem
using LineSearches
alg = Optim.BFGS(linesearch=LineSearches.BackTracking())
pythonplot()


# ------------------
# Load Kv function
# Computed in scripts/kv.jl, because it is used elsewhere too
@unpack K_shf = load(datadir("exp_pro", "simtra_kv.jld2"))
