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
    Sfac ~ truncated(Uniform(-15, 20))
    Cfac ~ truncated(Uniform(-9, 8))
    Rfac ~ truncated(Uniform(-4, 1))

    newpo = @set param.N = ConstPhysProp(exp(Nfac)*param.N(0))
    @reset newpo.S = ConstPhysProp(exp(Sfac)*param.S(0))
    @reset newpo.C = exp(Cfac)*param.C
    @reset newpo.Rs = exp(Rfac)*u"cm^2*Torr*hr/g"
    sol = solve(ODEProblem(newpo, tspan=400.0), Rodas3(), saveat=fitdat.t_Tf)

    sqerr = obj_expT(sol, fitdat, tweight=0)
    if !isnan(sqerr)
        Turing.@addlogprob! -sqerr
    else
        # @info "check" newpo.N.val newpo.S.val newpo.C.val newpo.Rs.val sol.t[end] sol[1,end]
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
    
alg = NUTS(10, 0.65; init_ϵ=1.0)
smallsamp = @time sample(case1_ns, alg, MCMCSerial(), 20, 1)
@profview sample(case1_ns, alg, MCMCSerial(), 20, 1)
plot(smallsamp)

alg = NUTS(50, 0.65)
samp = @time sample(case1_ns, alg, MCMCThreads(), 300, 12)

safesave(datadir("exp_pro", "turing_NS.jld2"), Dict("samp" => samp, "smallsamp" => smallsamp))

plot(samp)

Ns = exp.(reshape(samp[:Nfac], :))*basefit1.N(0)
Ss = exp.(reshape(samp[:Sfac], :))*basefit1.S(0)
prob = reshape(samp[:log_density], :)
scatter(Ns, Ss, marker_z=prob, scale=:log10, colorbar_title="log density")
plot!(xlabel="N", ylabel="S", )

