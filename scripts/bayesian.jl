using DrWatson
@quickactivate :FulbrightPores

using Turing

# -------------------
## Check methodology with an example from documentation
# using LinearAlgebra: I

# # Define Lotka-Volterra model.
# function lotka_volterra(du, u, p, t)
#     # Model parameters.
#     α, β, γ, δ = p
#     # Current state.
#     x, y = u
#     # Evaluate differential equations.
#     du[1] = (α - β * y) * x # prey
#     du[2] = (δ * x - γ) * y # predator
#     return nothing
# end

# # Define initial-value problem.
# u0 = [1.0, 1.0]
# p = [1.5, 1.0, 3.0, 1.0]
# tspan = (0.0, 10.0)
# prob = ODEProblem(lotka_volterra, u0, tspan, p)

# sol = solve(prob, Rosenbrock23(); saveat=0.1)
# odedata = Array(sol) + 0.8 * randn(size(Array(sol)))

# # Plot simulation and noisy observations.
# plot(sol; alpha=0.3)
# scatter!(sol.t, odedata'; color=[1 2], label="")

# @model function fitlv(data, prob)
#     # Prior distributions.
#     σ ~ InverseGamma(2, 3)
#     α ~ truncated(Normal(1.5, 0.5); lower=0.5, upper=2.5)
#     β ~ truncated(Normal(1.2, 0.5); lower=0, upper=2)
#     γ ~ truncated(Normal(3.0, 0.5); lower=1, upper=4)
#     δ ~ truncated(Normal(1.0, 0.5); lower=0, upper=2)

#     # Simulate Lotka-Volterra model. 
#     p = [α, β, γ, δ]
#     predicted = solve(prob, Rosenbrock23(); p=p, saveat=0.1)

#     # Observations.
#     for i in 1:length(predicted)
#         data[:, i] ~ MvNormal(predicted[i], σ^2 * I)
#     end

#     return nothing
# end
# @model function fitlv2(data, prob)
#     # Prior distributions.
#     σ ~ InverseGamma(2, 3)
#     α ~ truncated(Normal(1.5, 0.5); lower=0.5, upper=2.5)
#     β ~ truncated(Normal(1.2, 0.5); lower=0, upper=2)
#     γ ~ truncated(Normal(3.0, 0.5); lower=1, upper=4)
#     δ ~ truncated(Normal(1.0, 0.5); lower=0, upper=2)

#     # Simulate Lotka-Volterra model. 
#     p = [α, β, γ, δ]
#     predicted = solve(prob, Rosenbrock23(); p=p, saveat=0.1)

#     # Observations.
#     # for i in 1:length(predicted)
#     #     data[:, i] ~ MvNormal(predicted[i], σ^2 * I)
#     # end
#     sqerr = sum([sum(abs2, data[:,i] - predicted[i]) for i in eachindex(predicted)]) 
#     Turing.@addlogprob! -sqerr

#     return nothing
# end


# model = fitlv(odedata, prob)
# model2 = fitlv2(odedata, prob)

# # Sample 3 independent chains with forward-mode automatic differentiation (the default).
# testchain = @time sample(model, NUTS(), MCMCSerial(), 200, 3; progress=false)
# testchain2 = @time sample(model2, NUTS(), MCMCSerial(), 200, 3; progress=false)

# plot(testchain, c=:blue)
# plot!(testchain2, c=:red)

# ----------------
## Load fit results for my cases
fitdats = load(datadir("exp_pro", "simtra_sucrose_fitdata.jld2"))["fitdats"]
fits = load(datadir("exp_pro", "simtra_sucrose_fitresults.jld2"))
mcfits = fits["mc_fits"]

fitdat1 = fitdats[1]

basefit1 = mcfits[(Sigmoid_k, 1, (1,2,8))]

fit_guess = (2e6u"cm^-3", 5e-7u"cm^3", 2e-1u"m^0.5", 0.1u"cm^2*Torr*hr/g")


# -------------------------
# Turing
# https://turinglang.org/docs/tutorials/10-bayesian-differential-equations/

@model function porous_likelihood(fitdat, param)

    # σT ~ InverseGamma
    # Nfac ~ truncated(Normal(0.0, 2.0); lower=-8, upper=8)
    # Sfac ~ truncated(Normal(0.0, 2.0); lower=-8, upper=8)
    Nfac ~ truncated(Uniform(-6, 4))
    Sfac ~ truncated(Uniform(-12, 0))
    Cfac ~ truncated(Uniform(-4, 4))
    Rfac ~ truncated(Uniform(-2, 0))

    newpo = @set param.N = ConstPhysProp(exp10(Nfac)*fit_guess[1])
    @reset newpo.S = ConstPhysProp(exp10(Sfac)*fit_guess[2])
    @reset newpo.C = exp10(Cfac)*fit_guess[3]
    @reset newpo.Rs = exp10(Rfac)*u"cm^2*Torr*hr/g"
    sol = solve(ODEProblem(newpo, tspan=400.0), Rodas3())

    magfac = 10 # Make the squared error more important than priors

    sqerr = obj_expT(sol, fitdat, tweight=0) * magfac
    if !isnan(sqerr)
        Turing.@addlogprob! -sqerr
    else
        # @info "check" newpo.N.val newpo.S.val newpo.C.val newpo.Rs.val sol.t[end] sol[1,end]
        Turing.@addlogprob! -Inf
    end
    # fitdat.t_end ~ MvNormal(sol.t[end], 
    t_relerr = (sol.t[end]*u"hr" - fitdat.t_end)/fitdat.t_end
    t_sqerr = t_relerr^2 * magfac
    Turing.@addlogprob! t_sqerr
    # if sqerr
    # @info "see" sqerr t_sqerr
end

case1_params = porous_likelihood(fitdat1, basefit1)
    
# smallalg = NUTS(30, 0.65)
smallsamp = @time sample(case1_params, NUTS(), MCMCSerial(), 20, 1)
# @profview sample(case1_ns, alg, MCMCSerial(), 20, 1)
plot(smallsamp)

alg = NUTS(50, 0.65)
samp = @time sample(case1_params, alg, MCMCThreads(), 300, 4)

safesave(datadir("exp_pro", "case1_sampling.jld2"), Dict("samp" => samp, "smallsamp" => smallsamp))

# samp = load(datadir("exp_pro", "case1_sampling.jld2"))["samp"]

plot(samp)

# Ns = exp10.(reshape(samp[:Nfac], :))*basefit1.N(0)
# Ss = exp10.(reshape(samp[:Sfac], :))*basefit1.S(0)
# prob = reshape(samp[:log_density], :)

table = Table(map(eachindex(samp[:log_density])) do i
    (N = exp10(samp[:Nfac][i])*fit_guess[1],
     S = exp10(samp[:Sfac][i])*fit_guess[2],
     C = exp10(samp[:Cfac][i])*fit_guess[3],
     Rs = exp10(samp[:Rfac][i])*u"cm^2*Torr*hr/g",
     logprob = samp[:log_density][i],
     prob = exp(samp[:log_density][i]),
    )
end )
table = sort(table, by=x->x.logprob)

# Corner-plot style logic?
vars = [:N, :S, :C, :Rs]
for i in eachindex(vars)
    for j in i+1:length(vars)
    end
end

begin
pl1 = plot(u"mL^-1", u"mL", xlabel="N", ylabel="S")
@df table scatter!(:N, :S, marker_z=:logprob, scale=:log10, colorbar_title="log density")
end

# cutoff = [length(table.logprob)÷20] # 95%
HPD = table[length(table)÷20 : end] # 95% HPD, since table is sorted

#-----------------
# Trouble with a convex hull:
# you have to define what space it is in.
# If I use a log-log space, then the convex hull is different and I will have to take logs
# In the absence of any better justification, I'll just plan to provide#
# the scatter plot of the chain results, since that conveys the point pretty well


# using QHull

# nondim = hcat(ustrip.([HPD.N, HPD.S, HPD.C, HPD.Rs])...)
# hullNS = chull(nondim[:,1:2])
# polygon(nondim[hullNS.vertices,1], nondim[hullNS.vertices,2], )
# scatter!(nondim[:,1], nondim[:,2])
# # polyhed = polyhedron(vrep(nondim), QHull.Library())
# # removevredundancy!(polyhed)

# # plot(polyhed)
# # polyhed = polyhedron(hull)
# p = polyhedron(v, QHull.Library())
# # Removing redundant points, i.e. points which are in the interior of the convex hull

# plot(hull.simplices')

# marginalhist(Ns, Ss, xscale=:log10)
