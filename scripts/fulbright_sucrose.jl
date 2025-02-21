const OP = OptimizationProblem
using LineSearches
optalg = BFGS(linesearch=LineSearches.BackTracking())
default(:linewidth, 3)
default(:markersize, 5)

begin
Kv = load(datadir("exp_pro", "gravimetric_kv.jld2"))["Kv"]

dat1raw = CSV.read(datadir("exp_raw", "Fulbright2025", "2025-01-20-10_MFD_ISW.csv"), Table; header=8)
dat1pd = filter(x->x.Phase==4, dat1raw)
dat1 = map(x->microfd_fulbright_columnrename(x, dat1pd.CycleTime[1]), dat1pd)
for i in eachindex(dat1)[2:end]
    if dat1.t[i] < dat1.t[i-1]
        dat1.t[i:end] .+= 24u"hr"
    end
end
# plot(dat1raw.Phase)

pch1_pir_der2 = savitzky_golay(dat1.pch_pir, 91, 3, deriv=2, rate=1/u"minute").y
t_end1 = dat1.t[argmax(pch1_pir_der2[1u"hr" .< dat1.t .< 50u"hr" ])] + 1u"hr"

pch1 = RampedVariable(100u"mTorr")
# @df dat1 plot(:t, :pch_cm, ylim=(0, 200))
# @df dat1 plot!(:t, :pch_pir, ylim=(0, 200))
# plot!(pch1)
# tendplot!(t_end1)

Tsh1 = RampedVariable(uconvert.(u"K", [-47, -20]u"°C"), 0.5u"K/minute")
# @df dat1 plot(:t, :Tsh_i, xlim=(0, 2))
# plot!(Tsh1)

@df dat1 exptfplot(:t, :T1, :T2, :T3, :T4)

fitdat1 = @df dat1 PrimaryDryFit(:t, (
    :T1[:t .< 5u"hr"],
    :T2[:t .< 5u"hr"],
    :T3[:t .< 5u"hr"],
    :T4[:t .< 5.9u"hr"],
    ),
    t_end1
)
# plot(fitdat1)

dat2raw = CSV.read(datadir("exp_raw", "Fulbright2025", "2025-01-28-10_MFD_ISW.csv"), Table; header=8)
dat2pd = filter(x->x.Phase==4, dat2raw)
dat2 = map(x->microfd_fulbright_columnrename(x, dat2pd.CycleTime[1]), dat2pd)
for i in eachindex(dat2)[2:end]
    if dat2.t[i] < dat2.t[i-1]
        dat2.t[i:end] .+= 24u"hr"
    end
end

pch2_pir_der2 = savitzky_golay(dat2.pch_pir, 91, 3, deriv=2, rate=1/u"minute").y
t_end2 = dat2.t[argmax(pch2_pir_der2[1u"hr" .< dat2.t .< 50u"hr" ])] + 1u"hr"

pch2 = RampedVariable(100u"mTorr")
# @df dat2 plot(:t, :pch_cm, ylim=(0, 200))
# @df dat2 plot!(:t, :pch_pir, ylim=(0, 200))
# plot!(pch2)
# tendplot!(t_end2)

Tsh2 = RampedVariable(uconvert.(u"K", [-47, -10]u"°C"), 0.5u"K/minute")
# @df dat2 plot(:t, :Tsh_i, xlim=(0, 2))
# plot!(Tsh2)

# @df dat2 exptfplot(:t, :T1, :T2, :T3, :T4)

fitdat2 = @df dat2 PrimaryDryFit(:t, (
    :T1[:t .< 7u"hr"],
    :T2[:t .< 7u"hr"],
    :T3[:t .< 5u"hr"],
    :T4[:t .< 3u"hr"],
    ),
    t_end2
)
# plot(fitdat2)

dat3raw = CSV.read(datadir("exp_raw", "Fulbright2025", "2025-02-03-10_MFD_ISW.csv"), Table; header=8)
dat3pd = filter(x->x.Phase==4, dat3raw)
dat3 = map(x->microfd_fulbright_columnrename(x, dat3pd.CycleTime[1]), dat3pd)
for i in eachindex(dat3)[2:end]
    if dat3.t[i] < dat3.t[i-1]
        dat3.t[i:end] .+= 24u"hr"
    end
end

pch3_pir_der2 = savitzky_golay(dat3.pch_pir, 91, 3, deriv=2, rate=1/u"minute").y
t_end3 = dat3.t[argmax(pch3_pir_der2[1u"hr" .< dat3.t .< 50u"hr" ])] + 1u"hr"

pch3 = RampedVariable(100u"mTorr")
# @df dat3 plot(:t, :pch_cm, ylim=(0, 200))
# @df dat3 plot!(:t, :pch_pir, ylim=(0, 200))
# plot!(pch3)
# tendplot!(t_end3)

Tsh3 = RampedVariable(uconvert.(u"K", [-47, 0]u"°C"), 0.5u"K/minute")
# @df dat3 plot(:t, :Tsh_i, xlim=(0, 2))
# plot!(Tsh3)

# @df dat3 exptfplot(:t, :T1, :T2, :T3, :T4)

fitdat3 = @df dat3 PrimaryDryFit(:t, (
    :T1[:t .< 3u"hr"],
    :T2[:t .< 4u"hr"],
    :T3[:t .< 3.5u"hr"],
    :T4[:t .< 3.5u"hr"],
    ),
    t_end3
)
# plot(fitdat3)
end

blankplot_hrC()
plot!(fitdat1, c=1)
plot!(fitdat2, c=2)
plot!(fitdat3, c=3)




# ---------------
# blankplot_hrC(xlim=(0, 10), ylim=(-40, -30))
# @df dat1 exptfplot!(:t, :T1, :T2, :T3, :T4, c=1, labsuffix=", -20°C")
# @df dat2 exptfplot!(:t, :T1, :T2, :T3, :T4, c=2, labsuffix=", -10°C")
# savefig(plotsdir("sucrose_2cases.png"))

begin
rad_i, rad_o = get_vial_radii("10R")
A_p = π*rad_i^2
A_v = π*rad_o^2
fillvol = 2u"mL"
hf0 = fillvol/A_p |> u"cm"

Rp_g = RpFormFit(0.8u"cm^2*Torr*hr/g", 14u"cm*Torr*hr/g", 1u"1/cm")

c_solid = 0.05u"g/mL"
ρ_solution = 1u"g/mL"

d_typical = 50u"μm"
V_typical = 4/3*π*d_typical^3
N_typical = 1/V_typical
# Dummy values
N = ConstPhysProp(N_typical)
S = ConstPhysProp(V_typical) # = V_typ^2*N_typ
C = 2e-1u"m^0.5"
Rs = 0.5u"cm^2*Torr*hr/g"
fit_scl = (N(0), S(0), C, Rs)
k0 = 1e-14u"1/s"

pop1 = ParamObjPikal([
    (Rp_g, hf0, c_solid, ρ_solution),
    (Kv, A_v, A_p),
    (pch1, Tsh1),
])
pomc1 = ParamObjMC([
    (hf0, c_solid, ρ_solution),
    (Kv, A_v, A_p),
    (pch1, Tsh1),
    (ConstPhysProp(k0), C, N, S, Rs)
])

pop2 = @set pop1.Tsh = Tsh2
pop3 = @set pop1.Tsh = Tsh3
pomc2 = @set pomc1.Tsh = Tsh2
pomc3 = @set pomc1.Tsh = Tsh3
end

cases = Table([ (name="1", pop=pop1, pomc=pomc1, fitdat=fitdat1, dat=dat1),
                (name="2", pop=pop2, pomc=pomc2, fitdat=fitdat2, dat=dat2),
                (name="3", pop=pop3, pomc=pomc3, fitdat=fitdat3, dat=dat3),])

obj_Rpf = OptimizationFunction(obj_Rp, AutoForwardDiff())
@time obj_Rpf(fill(0.0, 3), (pop1, fitdat1))

obj_npf = OptimizationFunction(objn_Rp, AutoForwardDiff())


single_Rps = map(cases) do case
    (;name, pop, fitdat, dat) = case
    opt = solve(OP(obj_Rpf, [0.0, 0, -5], (pop, fitdat)), optalg, g_tol=1e-3)
    sol = gen_sol_Rp(opt.u, pop)
    @info "analyis" name opt sol
    Rp = sol.prob.p.Rp
    return (opt=opt, sol_p=sol, Rp=Rp)
end


cases_Rps = Table(cases, single_Rps)

pls = map(cases_Rps) do case 
    plT = blankplot_hrC()
    plot!(case.fitdat)
    modconvtplot!(case.sol_p)
    return plT
end
plot(pls..., layout=(3,1), size=(500, 1000))

function make_mcpass(casenums)
    subdat = filter(x->parse(Int, x.name)∈casenums, cases)
    (subdat.pomc, subdat.fitdat)
end

pass123 = make_mcpass([1, 2, 3])

myega_scale = (1e10u"s^-1", 1.6e3u"K", 1e2u"K")
k_f_myega(vec) = MYEGA_k((exp.(vec).*myega_scale))
objn_myega(fitprm, pf) = objn_mc(fitprm, (pf..., k_f_myega), tweight=1e-1, fit_scl=fit_scl)
objnf_myega = OptimizationFunction(objn_myega, AutoForwardDiff())
lb_myega = vcat(fill(-10.0, 4), fill(-100, 3))
ub_myega = vcat(fill(10.0, 4), fill( 100, 3))

sigm_scale = (1e-15u"s^-1", 1e-9u"s^-1", 0.8u"K", 239u"K")
function k_f_sigm(vec)
    p123 = exp.(vec[1:3] .-1).*sigm_scale[1:3]
    p4 = vec[4]*u"K" + sigm_scale[4]
    Sigmoid_k(p123..., p4)
end
objn_sigm(fitprm, pf) = objn_mc(fitprm, (pf..., k_f_sigm), tweight=1e-1, fit_scl=fit_scl)
objnf_sigm = OptimizationFunction(objn_sigm, AutoForwardDiff())
lb_sigm = vcat(fill(-10.0, 4), fill(-30, 4))
# lb_sigm[8] = -30
ub_sigm = vcat(fill(10.0, 4), fill( 30, 4))
# ub_sigm[8] = 30

@time opt123sigm = solve(OP(objnf_sigm, fill(0.1, 8), pass123, lb=lb_sigm, ub=ub_sigm), optalg, maxiters=50, time_limit=15)
p0_myega = [5.0, 1, 0, 0.1, 10, 1, 1]
@time opt123myega = solve(OP(objnf_myega, p0_myega, pass123, lb=lb_myega, ub=ub_myega), optalg, maxiters=50, time_limit=15, show_trace=true)

fits = Table((opt=[opt123sigm, opt123myega], 
              k_type=[k_f_sigm, k_f_myega],
              cases=[(1,2,3), (1,2,3)]))
fits = Table(fits, (kf= map(x->x.k_type(x.opt[5:end]), fits),))

function gen_pbmsols(opt, k_f_gen)
    map(cases) do case
        sol = gen_sol_mc_arb(opt.u, case.pomc, k_f_gen; fit_scl=fit_scl)
        return (sol_mc=sol,k_f=sol.prob.p.k_f)
    end
end

cases_both = Table(cases, single_Rps, gen_pbmsols(opt123myega, k_f_myega))
# cases_both = Table(cases, single_Rps, gen_pbmsols(opt123sigm, k_f_sigm))

begin
pls_T = map(cases_both) do case
    pl = blankplot_hrC(legend=:outerright, title=case.name)
    plot!(case.fitdat)
    modconvtplot!(case.sol_mc, case.sol_p, linealpha=0.7, markeralpha=0, labels=["PBM: M 1&2&3" "3PP $(case.name)"])
    plot!(xlim=(0, 16), ylim=(-50, -28))
    return pl
end
# plot!(pls_T[5], legend=:bottomright)
plot(pls_T..., layout=(3,1), size=(500, 600))
end
savefig(plotsdir("fulbright_sucrose_T.svg"))
savefig(plotsdir("fulbright_sucrose_T.pdf"))

l = range(0u"cm", 0.5u"cm", length=100)
begin
poreevolplot(cases_both[1].sol_mc, label="PBM 1, 1&2&3", marker=:o)
poreevolplot!(cases_both[2].sol_mc, label="PBM 2, 1&2&3", marker=:square)
poreevolplot!(cases_both[3].sol_mc, label="PBM 3, 1&2&3", marker=:d)
plot!(subplot=4, l, cases_both.Rp[1].(l), ls=:dash, c=1, label="3PP 1")
plot!(subplot=4, l, cases_both.Rp[2].(l), ls=:dash, c=2, label="3PP 2")
plot!(subplot=4, l, cases_both.Rp[3].(l), ls=:dash, c=3, label="3PP 3")
end

res_fit = cases_both[1].sol_mc.prob.p
@info "Best fit params:" res_fit.N res_fit.S res_fit.C res_fit.Rs

function eval_about_myega(new, opt, pass)
    objn_myega(new .* opt.u, pass)
end
function eval_many(ranges...)
    ndims = length(ranges)
    eval_points = [[u..., fill(1.0, 7-ndims)...] for u in Iterators.product(ranges...)]
    eval_about_myega.(eval_points, [opt123myega], [pass123])
end

# spread = [[exp(i), exp(j), fill(1.0,5)...] for i in -5:5, j in -5:5]

integ_dims(dims...) = dropdims(sum(probs, dims=dims), dims=dims)
# @time sqerrs = eval_about_myega.(spread, [opt123myega], [pass123])
begin
Nr = range(-10, 10, length=11)
Sr = range(-2, 4, length=12)
Cr = range(-1, 4, length=5)
Rsr = range(0, 3.5, length=4)
@time sqerrs = eval_many(Nr, Sr, Cr, Rsr)
minsqerr = minimum(filter(!isnan, sqerrs)) 
opt123myega.minimum
probs = exp.(-map(x->(isnan(x) ? Inf : x-minsqerr), sqerrs))
probs ./= sum(probs)
end
contourf(integ_dims(1, 2), cmap=:thermal, cscale=:log10, lw=0)
contourf(integ_dims(3, 4), cmap=:thermal, cscale=:log10, lw=0)
# contourf(sqerrs[:,:,1,1], cmap=:thermal, lw=0)



# ## Bayesian, a la DynamicHMC

# struct MYEGA_prob1{P, F}
#     po::P
#     fitdat::F
# end

# function (prob::MYEGA_prob1)(θ)
#     lim = [10, 3, 2, 2, 5, 5, 5]
#     if any(θ .< -lim) || any(θ .> lim)
#         return -Inf
#     end
#     rtype = eltype(θ)
#     sqerr::rtype = objn_myega(opt123myega.u .+ θ, (prob.po, prob.fitdat))
#     if isnan(sqerr)
#         return rtype(-Inf)
#     else
#         return -sqerr
#     end
# end

# bayes123 = MYEGA_prob1(pass123...)

# # logprob = LogDensityProblem(bayes123)
# using LogDensityProblems, LogDensityProblemsAD
# LogDensityProblems.dimension(::MYEGA_prob1) = 7
# LogDensityProblems.logdensity(prob::MYEGA_prob1, θ) = -prob(θ)
# LogDensityProblems.capabilities(bayes123)

# import Random
# using DynamicHMC
# ∇l = ADgradient(:ForwardDiff, bayes123)
# @time bayes123(fill(-3.0, 7))
# # @code_warntype LogDensityProblems.logdensity_and_gradient(∇l, fill(1.0, 7))
# # @time ∇l(fill(1.0, 7))
# # @time results = mcmc_with_warmup(Random.default_rng(), ∇l, 10; reporter = NoProgressReport())
# @time results = mcmc_with_warmup(Random.default_rng(), ∇l, 10, initialization=(q=fill(0.0, 7),ϵ=0.01))
# @time results = mcmc_with_warmup(Random.default_rng(), ∇l, 1000, initialization=(q=fill(0.0, 7),ϵ=0.01), warmup_stages=fixed_stepsize_warmup_stages())

# safesave(datadir("exp_pro", "suc_mcmc.jld2"), @strdict results sqerrs Nr Sr Cr Rsr probs)