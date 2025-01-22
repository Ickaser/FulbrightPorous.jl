using DrWatson
@quickactivate :FulbrightPores

using Turing

# ----------------
## Load fit results
fitdats = load(datadir("exp_pro", "simtra_sucrose_fitdata.jld2"))["fitdats"]
fits = load(datadir("exp_pro", "simtra_sucrose_fitresults.jld2"))
mcfits = fits["mc_fits"]

fitdat1 = fitdats[1]

basefit1 = mcfits[(Sigmoid_k, 1, (1,2,8))]



# -------------------------
# Turing
# https://turinglang.org/docs/tutorials/10-bayesian-differential-equations/

@model function turing_NS(fitdat, param)

    # σT ~ InverseGamma
    # Nfac ~ truncated(Normal(0.0, 2.0); lower=-8, upper=8)
    # Sfac ~ truncated(Normal(0.0, 2.0); lower=-8, upper=8)
    Nfac ~ truncated(Uniform(-8, 8))
    Sfac ~ truncated(Uniform(-8, 8))
    Cfac ~ truncated(Uniform(-8, 8))
    Rfac ~ truncated(Uniform(-8, 8))

    newpo = @set param.N = ConstPhysProp(exp(Nfac)*param.N(0))
    @reset newpo.S = ConstPhysProp(exp(Sfac)*param.S(0))
    @reset newpo.C = ConstPhysProp(exp(Cfac)*param.C)
    @reset newpo.Rs = ConstPhysProp(exp(Rfac)*param.Rs)
    sol = solve(ODEProblem(newpo), Rodas3())

    sqerr = obj_expT(sol, fitdat, tweight=0)
    if !isnan(sqerr)
        Turing.@addlogprob! -sqerr
    else
        Turing.@addlogprob! -Inf
    end
    # errs = err_expT(sol, fitdat)
    # ntot = length(errs)
    # skip = ntot ÷ 100 # Want roughly 100 data points to not overwhelm things
    # downsampled = errs[1:skip:end-1]
    # for i in eachindex(downsampled)
    #     downsampled[i] ~ Normal(0.0, 1.0)
    # end
    # errs[end] ~ Normal(0.0, 10.0)
end

case1_ns = turing_NS(fitdat1, basefit1)
    
@time sample(case1_ns, NUTS(), MCMCSerial(), 20, 1)

samp = @time sample(case1_ns, NUTS(), MCMCSerial(), 1000, 10)

plot(samp)

# Ns = exp.(reshape(samp[:Nfac], :))*fit_scl_porous[1]
# Ss = exp.(reshape(samp[:Sfac], :))*fit_scl_porous[2]
# prob = reshape(samp[:log_density], :)
# scatter(Ns, Ss, marker_z=prob, scale=:log10, )

# If Turing is annoying, could try DynamicHMC.jl.
# https://www.tamaspapp.eu/DynamicHMC.jl/dev/worked_example/