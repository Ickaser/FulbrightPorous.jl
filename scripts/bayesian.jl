using DrWatson
@quickactivate :FulbrightPores

using Turing

# ----------------
## Load fit results for my cases
fitdats = load(datadir("exp_pro", "simtra_sucrose_fitdata.jld2"))["fitdats"]
fits = load(datadir("exp_pro", "simtra_sucrose_fitresults.jld2"));
mcfits = fits["mc_fits"]

fitdat1 = fitdats[1]

basefit1 = mcfits[(MYEGA_k, 1, (1,2,8))]
pos128 = mcfits[(MYEGA_k, 1, (1,2,8))], mcfits[(MYEGA_k, 2, (1,2,8))], mcfits[(MYEGA_k, 8, (1,2,8))]
fits128 = fitdats[1], fitdats[2], fitdats[8]

# porous_scale = (2e6u"cm^-3", 5e-7u"cm^3", 2e-1u"m^0.5", 0.1u"cm^2*Torr*hr/g")
porous_scale = FulbrightPores.fit_scl_porous
basek = basefit1.k_f
myega_scale = (basek.k0, basek.K, basek.C)
# myega_scale = (3e6u"s^-1", 5e2u"K", 1e3u"K")

k_f_myega_nolog(vec) = MYEGA_k(vec.*myega_scale)

baseprob = ODEProblem(basefit1, tspan=(0.0, 1000.0))

# -------------------------

@model function myega_likelihood(param, fitdat)

    magfac = 10 # Make the squared error more important than priors

    # σT ~ InverseGamma
    # Nfac ~ truncated(Normal(0.0, 2.0); lower=-8, upper=8)
    # Sfac ~ truncated(Normal(0.0, 2.0); lower=-8, upper=8)
    Nfac ~ Uniform(-6, 4)
    Sfac ~ Uniform(-12, 0)
    Cfac ~ Uniform(-4, 4)
    Rfac ~ Uniform(-2, 1)
    kfac ~ Uniform(-2, 2)
    Kfac ~ Uniform(-2, 2)
    θfac ~ Uniform(-2, 2)

    fitprm = exp10.([Nfac, Sfac, Cfac, Rfac, kfac, Kfac, θfac])

    # k_args = exp10.([kfac, Kfac, θfac]) .* fit_guess[5:7]

    newpo = @set param.N = ConstPhysProp(fitprm[1]*porous_scale[1])
    @reset newpo.S = ConstPhysProp(fitprm[2]*porous_scale[2])
    @reset newpo.C = fitprm[3]*porous_scale[3]
    @reset newpo.Rs = fitprm[4]*porous_scale[4]
    @reset newpo.k_f = k_f_myega_nolog(fitprm[5:end])

    newprob = remake(baseprob, p=newpo)
    
    sol = solve(newprob, Rodas3(), saveat=ustrip.(u"hr", fitdat.t))

    sqerr = magfac * obj_expT(sol, fitdat, tweight=1.0) 
    # sqerr = magfac * objn_mc(fitprm, (params, fitdats, k_f_myega_nolog), tweight=1e-1) 

    if !isnan(sqerr)
        Turing.@addlogprob! -sqerr
    else
    #     # @info "check" newpo.N.val newpo.S.val newpo.C.val newpo.Rs.val sol.t[end] sol[1,end]
        Turing.@addlogprob! -Inf
    end
    # fitdat.t_end ~ MvNormal(sol.t[end], 
    # t_relerr = (sol.t[end]*u"hr" - fitdat.t_end)/fitdat.t_end
    # t_sqerr = t_relerr^2 * magfac
    # Turing.@addlogprob! t_sqerr
end

case1_params = myega_likelihood(basefit1, fitdat1)
# case128_params = myega_likelihood(pos128, fits128)
    
smallalg = NUTS(10, 0.65; init_ϵ=0.01)
smallsamp = @time sample(case1_params, smallalg, MCMCSerial(), 20, 1)
# smallsamp128 = @time sample(case128_params, smallalg, MCMCThreads(), 20, 4)
# @profview sample(case1_params, smallalg, MCMCSerial(), 20, 1)
plot(smallsamp)

alg = NUTS(75, 0.65; init_ϵ=0.01)
samp = @time sample(case1_params, alg, MCMCThreads(), 500, 4)
safesave(datadir("exp_pro", "case1_myega_sampling.jld2"), Dict("samp" => samp, "smallsamp" => smallsamp))
# samp128 = @time sample(case128_params, alg, MCMCThreads(), 300, 4)

# safesave(datadir("exp_pro", "case1_myega_sampling.jld2"), Dict("samp" => samp, "smallsamp" => smallsamp))
# safesave(datadir("exp_pro", "case128_myega_sampling.jld2"), @strdict smallsamp smallsamp128 samp )

# samp = load(datadir("exp_pro", "case1_sampling.jld2"))["samp"]

plot(samp)

# Ns = exp10.(reshape(samp[:Nfac], :))*basefit1.N(0)
# Ss = exp10.(reshape(samp[:Sfac], :))*basefit1.S(0)
# prob = reshape(samp[:log_density], :)

table = Table(map(eachindex(samp[:log_density])) do i
    (N = exp10(samp[:Nfac][i])*porous_scale[1],
     S = exp10(samp[:Sfac][i])*porous_scale[2],
     C = exp10(samp[:Cfac][i])*porous_scale[3],
     Rs = exp10(samp[:Rfac][i])*porous_scale[4],
     k0 = exp10(samp[:kfac][i])*myega_scale[1],
     K = exp10(samp[:Kfac][i])*myega_scale[2],
     θ = exp10(samp[:θfac][i])*myega_scale[3],
     logprob = samp[:log_density][i],
     prob = exp(samp[:log_density][i]),
    )
end )
table = sort(table, by=x->x.logprob)

# Corner-plot style logic?
vars = [:N, :S, :C, :Rs, :k0, :K, :θ]
pl = plot(layout=(6,6))
for i in eachindex(vars)[begin:end-1]
    for j in i+1:length(vars)
        plot!(subplot=6*(i-1)+j-1, getproperty(table, vars[i]), getproperty(table, vars[j]))
    end
end

begin
# pl1 = plot(u"mL^-1", u"mL", xlabel="N", ylabel="S")
@df table scatter(:Rs, :k0, marker_z=:logprob, scale=:log10, colorbar_title="log density")
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
