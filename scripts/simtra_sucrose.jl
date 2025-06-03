using DrWatson
@quickactivate :FulbrightPores

const OP = OptimizationProblem
using LineSearches
optalg = BFGS(linesearch=LineSearches.BackTracking())
default(:linewidth, 3)
default(:markersize, 5)
using UnitfulLatexify

begin
# ------------------
# Load Kv function
# Computed in scripts/kv.jl, because it is used elsewhere too
Kshf = load(datadir("exp_pro", "simtra_kv.jld2"))["K_shf"]

# --------------------------
# Cycle 1
# Labeled as case 4 in Table 8
dat1raw = CSV.read(datadir("exp_raw", "Simtra_Rp", "Suc_Cycle1.csv"), Table)
dat1 = map(simtra_rename, dat1raw)

t_end1 = identify_pd_end_2der(dat1, 391, 10u"hr", 85u"hr")

pch1 = RampedVariable(100u"mTorr")
# @df dat1 plot(:t, :pch_cm, ylim=(0, 200))
# @df dat1 plot(:t, :pch_pir, ylim=(0, 200))
# plot!(pch1)
# tendplot!(t_end1)

# There's a slight delay before shelf temperature ramps
Tsh1 = RampedVariable(uconvert.(u"K", [-40, -40, -25]u"°C"), [1u"K/minute", 0.5u"K/minute"], [23u"minute"])
# @df dat1 plot(:t, :Tsh, xlim=(0, 2))
# plot!(Tsh1)

fitdat1 = PrimaryDryFit(dat1.t, (
    filter(x->(!isnan(x) && x<-32.5u"°C"), dat1.T1),
    filter(!isnan, dat1.T2[dat1.t .< 60u"hr"]),
    # collect(filter(!isnan, dat1.T3)),
    ),
    t_end1
)
# plot(fitdat1)

# ------------------------------
# Cycle 2

dat2raw = CSV.read(datadir("exp_raw", "Simtra_Rp", "Suc_Cycle2.csv"), Table)
dat2 = map(simtra_rename, dat2raw)
filter!(x->x.t .> 0u"hr", dat2)

# @df dat1 exptfplot(:t, :T1, :T2, :T3)
# @df dat2 exptfplot!(:t, :T1, :T2, :T3, c=:red)
# plot!(ylim=(-40, -30))

t_end2 = identify_pd_end_2der(dat2, 391, 10u"hr", 120u"hr")
pch2 = RampedVariable(100u"mTorr")
# @df dat2 plot(:t, :pch_cm, ylim=(0, 200))
# @df dat2 plot(:t, :pch_pir, ylim=(0, 200))
# plot!(pch2)
# tendplot!(t_end2)

fitdat2 = @df dat2 PrimaryDryFit(:t, (
    # collect(filter(!isnan, :T1)),
    filter(!isnan, :T2[:t .< 85u"hr"]),
    filter(!isnan, :T3),),
    t_end2)
# plot(fitdat2)

Tsh2 = Tsh1 - 5u"K"
# @df dat2 plot(:t, :Tsh, xlim=(0, 2))
# plot!(Tsh2)

# -------------------------------  
# Cycle 3
dat3raw = CSV.read(datadir("exp_raw", "Simtra_Rp", "Suc_Cycle3.csv"), Table)
dat3 = map(simtra_rename, dat3raw)

pch3 = RampedVariable(60u"mTorr")
t_end3 = identify_pd_end_2der(dat3, 391, 10u"hr", 68u"hr")
# @df dat3 plot(:t, :pch_cm, ylim=(0, 200))
# @df dat3 plot!(:t, :pch_pir, ylim=(0, 200))
# plot!(pch3)
# tendplot!(t_end3)

Tsh3 = RampedVariable(uconvert.(u"K", [-40, -25]u"°C"), 0.5u"K/minute")
# @df dat3 plot(:t, :Tsh, xlim=(0, 2))
# plot!(Tsh3)
fitdat3 = @df dat3 PrimaryDryFit(:t, (filter(!isnan, :T1), filter(!isnan, :T2), :T3[:t.<45u"hr"]), t_end3)
# plot(fitdat3)

# -----------------
# Cycle 5
dat5raw = CSV.read(datadir("exp_raw", "Simtra_Rp", "Suc_Cycle5.csv"), Table)
dat5 = map(simtra_rename, dat5raw)

pch5 = RampedVariable(50u"mTorr")
t_end5 = identify_pd_end_2der(dat5, 391, 10u"hr", 50u"hr")
# @df dat5 plot(:t, :pch_cm, ylim=(0, 200))
# @df dat5 plot!(:t, :pch_pir, ylim=(0, 200))
# plot!(pch5)
# tendplot!(t_end5)

Tsh5 = RampedVariable(uconvert.(u"K", [-40, -40, -20]u"°C"), [1u"K/minute", 0.5u"K/minute"], [12u"minute"])
# @df dat5 plot(:t, :Tsh_i, xlim=(0, 2))
# plot!(Tsh5)

fitdat5 = @df dat5 PrimaryDryFit(:t, (filter(!isnan, :T1), filter(!isnan, :T3), :T2[:t.<40u"hr"]), t_end5)
# plot(fitdat5)

# -----------------
# Cycle 8, high temperature
dat8raw = CSV.read(datadir("exp_raw", "Simtra_Rp", "Suc_Cycle8.csv"), Table)
dat8 = map(simtra_rename, dat8raw)
# @df dat8 exptfplot(:t, :T1, :T2, :T3)
# @df dat8 plot!(:t, :Tsh_i, xlim=(0, 20), ylim=(-40, -20))

pch8 = RampedVariable(100u"mTorr")
t_end8 = identify_pd_end_2der(dat8, 391, 10u"hr", 41u"hr")
# @df dat8 plot(:t, :pch_cm, ylim=(0, 200))
# @df dat8 plot(:t, :pch_pir, ylim=(0, 200))
# plot!(pch8)
# tendplot!(t_end8)

Tsh8 = RampedVariable(uconvert.(u"K", [-40, -40, -10]u"°C"), [1u"K/minute", 0.5u"K/minute"], [4u"minute"])
# @df dat8 plot(:t, :Tsh_i, xlim=(0, 2))
# plot!(Tsh8)

fitdat8 = @df dat8 PrimaryDryFit(:t, (
    filter(!isnan, :T1),
    filter(!isnan, :T2[:t .< 17u"hr"]),
    filter(!isnan, :T3[:t .< 17u"hr"]),),
    t_end8);
# plot(fitdat8)
end

@info "All data loaded"

# -------------
# Parameter setup

begin
rad_i = 1.1u"cm"
rad_o = 1.2u"cm"
Ap = π*rad_i^2
Av = π*rad_o^2
fillvol = 5u"mL"
hf0 = fillvol/Ap

Rp_g = RpFormFit(0.8u"cm^2*Torr*hr/g", 14u"cm*Torr*hr/g", 1u"1/cm")

csolid = 0.04u"g/mL"
ρsolution = 1u"g/mL"

d_typical = 50u"μm"
V_typical = 4/3*π*d_typical^3
N_typical = 1/V_typical
# Dummy values
N = ConstPhysProp(N_typical)
S = ConstPhysProp(V_typical) # = V_typ^2*N_typ
C = 2e-1u"m^0.5"
Rs = 0.5u"cm^2*Torr*hr/g"
k0 = 1e-14u"1/s"

pop1 = ParamObjPikal([
    (Rp_g, hf0, csolid, ρsolution),
    (Kshf, Av, Ap),
    (pch1, Tsh1),
])
pomc1 = ParamObjMC([
    (hf0, csolid, ρsolution),
    (Kshf, Av, Ap),
    (pch1, Tsh1),
    (ConstPhysProp(k0), C, N, S, Rs)
])

function new_Tp(po, Tsh, pch)
    po_n = deepcopy(po)
    @reset po_n.pch = pch
    @reset po_n.Tsh = Tsh
    return po_n
end
pop2 = new_Tp(pop1, Tsh2, pch2)
pop3 = new_Tp(pop1, Tsh3, pch3)
pop5 = new_Tp(pop1, Tsh5, pch5)
pop8 = new_Tp(pop1, Tsh8, pch8)
pomc2 = new_Tp(pomc1, Tsh2, pch2)
pomc3 = new_Tp(pomc1, Tsh3, pch3)
pomc5 = new_Tp(pomc1, Tsh5, pch5)
pomc8 = new_Tp(pomc1, Tsh8, pch8)
dat1 = Table(dat1, (ld=FulbrightPores.calc_l(dat1, pop1.hf0)))
dat2 = Table(dat2, (ld=FulbrightPores.calc_l(dat2, pop2.hf0)))
dat3 = Table(dat3, (ld=FulbrightPores.calc_l(dat3, pop3.hf0)))
dat5 = Table(dat5, (ld=FulbrightPores.calc_l(dat5, pop5.hf0)))
dat8 = Table(dat8, (ld=FulbrightPores.calc_l(dat8, pop8.hf0)))
end

cases = Table([ (name="1", pop=pop1, pomc=pomc1, fitdat=fitdat1, dat=dat1),
                (name="2", pop=pop2, pomc=pomc2, fitdat=fitdat2, dat=dat2),
                (name="3", pop=pop3, pomc=pomc3, fitdat=fitdat3, dat=dat3),
                (name="5", pop=pop5, pomc=pomc5, fitdat=fitdat5, dat=dat5),
                (name="8", pop=pop8, pomc=pomc8, fitdat=fitdat8, dat=dat8),
                ])

trans_Rp = Rp_transform_basic(0.8u"cm^2*Torr*hr/g", 14.0u"cm*Torr*hr/g", 1e0u"1/cm")
obj_Rpf = OptimizationFunction((x,y)->obj_pd(x,y,tweight=1e-2), AutoForwardDiff(chunksize=3))
# obj_Rpf = OptimizationFunction(obj_pd, AutoForwardDiff(chunksize=3))
@time obj_Rpf([0.0, 1.0, 0.0], (trans_Rp, pop1, fitdat1))

obj_npf = OptimizationFunction(objn_pd, AutoForwardDiff(chunksize=3))

single_Rps = map(cases) do case
    (;name, pop, fitdat, dat) = case
    u0 = [0., 0., 0.]
    @time opt = solve(OP(obj_Rpf, u0, (trans_Rp, pop, fitdat)), optalg)#, g_tol=1e-6)
    Rp = RpFormFit(transform(trans_Rp, opt.u).Rp...)
    sol = gen_sol_pd(opt.u, trans_Rp, pop, fitdat)
    @info "analyis" name opt Rp
    hRps = map(x->calc_hRp_T(case.pop, case.fitdat, i=x), 1:length(case.fitdat.Tfs))
    return (opt=opt, sol_p=sol, Rp, hRps)
end;

cases_Rps = Table(cases, single_Rps)

pls = map(cases_Rps) do case 
    plT = blankplot_hrC()
    plot!(case.fitdat, title=case.name)
    modconvtplot!(case.sol_p)
    plot!(legend = :bottomright)
    # plot!(case.pop.Tsh, tmax=case.fitdat.t_end, c=:black, label=L"T_\mathrm{sh}")
    return plT
end
plot(pls..., layout=(3,2), size=(800, 1000))

l = range(0u"cm", hf0, length=100)
pls = map(cases_Rps) do case
    plr = blankplot_lRp()
    plot!(case.dat.ld, case.dat.Rp, label="TDLAS")
    for i in 1:length(case.fitdat.Tfs)
        hdT, RpT = case.hRps[i]
        plot!(hdT, RpT, c=i+1, label="Rp from T, $i")
    end
    plot!(l, case.Rp.(l), c=5, label="3PP fit directly to T")
    plot!(ylim=(0, 10), legend = :bottomright)
    return plr
end
plot(pls..., layout=(3,2), size=(800, 1000))

function make_Rp_pass(casenums)
    subdat = filter(x->parse(Int, x.name)∈casenums, cases_Rps)
    (subdat.pop, subdat.fitdat)
end
function make_mc_pass(casenums)
    subdat = filter(x->parse(Int, x.name)∈casenums, cases)
    (subdat.pomc, subdat.fitdat)
end

rpfit_cases = Table(cases = [(1,), (2,), (3,), (5,), (8,), (1,8), (1,2,8)],
        p0 = fill(zeros(3), 7),
        trans = fill(trans_Rp, 7))
rpfits = map(rpfit_cases) do case
    @time solve(OP(obj_npf, case.p0, (case.trans, make_Rp_pass(case.cases)...)), optalg;)
end

# -----------------------------
# Fit k by interpolating at four values
# interp_nt = (kl1 = TV.Identity(),
#             kl2 = TV.Identity(),
#             kl3 = TV.Identity(),
#             kl4 = TV.Identity())
const TV = TransformVariables
interp_nt = (kl1 = TV.Identity(),
            kl2 = TVExp(),
            kl3 = TVExp(),
            kl4 = TVExp())
function FulbrightPores.k_f_nt(nt::NamedTuple{(:kl1, :kl2, :kl3, :kl4)})
    # return InterpK([nt...], [-45.0, -38, -32, -28]u"°C" .|> u"K", 1e-14u"s^-1")
    # if nt.kl1 < nt.kl2 || nt.kl2 < nt.kl3 || nt.kl3 < nt.kl4
    #     return 
    # end
    return InterpK(cumsum([nt...]), [-45.0, -38, -32, -22]u"°C" .|> u"K", 1e-14u"s^-1")
end
const TV = TransformVariables
porous_nt = (N = ConstWrapTV() ∘ TVScale(3e7u"cm^-3") ∘ TVExp(),
               S = ConstWrapTV() ∘ TVScale(2e-5u"cm^3") ∘ TVExp(), # If I could enforce S > ϕ/N, would be great...
               C = TVScale(4e-1u"m^0.5") ∘ TVExp(),
               Rs = TVScale(0.5u"cm^2*Torr*hr/g") ∘ TVShift(0.01) ∘ TVLogistic()
               )
myega_nt = (k0 = TVScale(1e11u"s^-1") ∘ TVExp(),
            K1 = TVScale(5e4u"K") ∘ TVShift(0.001) ∘ TVLogistic() ∘ TVShift(TV.logit(0.05)),
            K2 = TVScale(6e2u"K") ∘ TVShift(0.01) ∘ TVLogistic() ,)
# plot(Tran, FulbrightPores.k_f_nt(transform(as(myega_nt), zeros(3))).(Tran), yscale=:log10, label="MYEGA, def")
sigma_nt = (k0 = TVScale(1e-15u"s^-1") ∘ TVExp(),
            kf = TVScale(1e-8u"s^-1") ∘ TVExp(),
            K1 = TVScale(0.8u"K") ∘ TVExp(),
            K2 = TVScale(1.0u"K") ∘ TVShift(216.0) ∘ TVScale(45.0) ∘ TVLogistic() ∘TVScale(0.1) ,)
trans_interp = as(merge(porous_nt, (;k_f=as(interp_nt))))
trans_myega = as(merge(porous_nt, (;k_f=as(myega_nt))))
trans_sigma = as(merge(porous_nt, (;k_f=as(sigma_nt))))

obj_nmc = OptimizationFunction((x,y)->objn_pd(x,y, tweight=1e-2), AutoForwardDiff(chunksize=4))


interp_guess = [-1.0, 1, 0, 0, -1, 1.5, 0, 2]
# plot(u"°C", u"s^-1")
# kineticplot!(transform(trans_interp, interp_guess),  label="Interp guess")
myega_guess = [0.1, 0.1, 0, 0, -10, 0, 0]
# kineticplot(transform(trans_myega, myega_guess), label="MYEGA guess")
sigma_guess = [-1.0, 1, 0, 0, 0, 0, 0, 0]
# kineticplot(transform(trans_sigma, sigma_guess), label="MYEGA guess")

mcfit_cases = Table(cases = repeat([(1,), (8,), (1, 8), (1, 2, 8)], 3),
        p0 = vcat(fill(interp_guess, 4), fill(myega_guess, 4), fill(sigma_guess, 4)),
        trans = vcat(fill(trans_interp, 4), fill(trans_myega, 4), fill(trans_sigma, 4)))

mcfits = map(mcfit_cases) do case
    @time solve(OP(obj_nmc, case.p0, (case.trans, make_mc_pass(case.cases)...)), optalg;)
    # @time solve(OP(obj_nmc, case.p0, (case.trans, make_mc_pass(case.cases)...)), optalg; show_trace=true)
end

begin
plot(u"°C", u"s^-1")
for i in [3,4,7,8,11,12]
    k = k_f_nt(transform(mcfit_cases[i].trans, mcfits[i].u).k_f)
    kineticplot!(k, label="$(nameof(typeof(k))) fit to $(mcfit_cases[i].cases)")
end
# kineticplot!(transform(trans_interp, interp_guess),  label="Interp guess")
plot!(ylim=(1e-20, 1e-0), legend=:bottomright)
end


function get_sol_mc(casenum, fitnum)
    case = filter(x->parse(Int, x.name) == casenum, cases_Rps)[1]
    sol = gen_sol_pd(mcfits[fitnum].u, mcfit_cases[fitnum].trans, case.pomc)
end
function get_sol_Rp(casenum, fitnum)
    case = filter(x->parse(Int, x.name) == casenum, cases_Rps)[1]
    sol = gen_sol_pd(rpfits[fitnum].u, rpfit_cases[fitnum].trans, case.pop)
end

# Choose a favored fit for plotting, etc.
# fitnum = 7 # 1&8, MYEGA
fitnum = 11 # 1&8, sigmoid

begin
pls = map(cases_Rps, [1, 2, 3, 5, 8]) do case, i
    pl = blankplot_hrC(title=i, legend=:bottomright)
    plot!(case.fitdat)
    modconvtplot!(get_sol_mc(i, fitnum), label="PBM, fit to $(mcfit_cases[fitnum].cases)")
    modconvtplot!(case.sol_p, c=:green, alpha=0.6, label="3PP, fit to $(case.name)")
    modconvtplot!(get_sol_Rp(i, 6), c=:purple, label="3PP, fit to $(rpfit_cases[6].cases)")
    # plot!(xlim=(-1, 120), ylim=(-45, -28)) 
    return pl
end
plot(pls..., layout=(3,2), size=(800, 1000))
end
savefig(plotsdir("simtra_suc_fitval.svg"))
savefig(plotsdir("simtra_suc_fitval.pdf"))

safesave(datadir("exp_pro", "simtra_sucrose_fitdata.jld2"), Dict("fitdats"=>cases.fitdat))
safesave(datadir("exp_pro", "simtra_sucrose_fitresults.jld2"), @strdict mcfit_cases mcfits rpfit_cases rpfits)

begin
plot!(palette=:default)
rpattr = (seriestype=:samplemarkers, linestyle=:dash)
poreevolplot( get_sol_mc(1, fitnum), label="PBM 1&2&8, Cycle 1", c=1, marker=:circle)
poreevolplot!(get_sol_mc(2, fitnum), label="PBM 1&2&8, Cycle 2", c=2, marker=:square)
poreevolplot!(get_sol_mc(8, fitnum), label="PBM 1&2&8, Cycle 8", c=3, marker=:utriangle)
plot!(l, cases_Rps[1].Rp.(l), subplot=4, label="3PP 1"; c=4, marker=:circle, rpattr...)
plot!(l, cases_Rps[2].Rp.(l), subplot=4, label="3PP 2"; c=4, marker=:square, rpattr...)
plot!(l, cases_Rps[5].Rp.(l), subplot=4, label="3PP 8"; c=4, marker=:utriangle, rpattr...)
    
end

# -----------------
# Inspect fit results in Rp space rather than temperature space

begin
pls = map(cases_Rps[[1, 2, 4, 5]], [1,2,5,8]) do case, i
    pl = blankplot_lRp(ylims=(0, 12), xlims=(0, 1.3), legend=:bottomright)
    plot!(title="Cycle $i, \$T_{sh}=$(float(case.pop.Tsh(Inf*u"s")|>u"°C"))\$")
    pbmtdlas!(get_sol_mc(i, fitnum), case.dat, case.Rp)
    i ∈ [1, 2] && plot!(xlabel="", xformatter=:none)
    i ∈ [2, 8] && plot!(ylabel="", yformatter=:none)
    pl
end
plot(pls..., layout=(2,2), size=(600,500))
end
savefig(plotsdir("suc18_Rp_val.svg"))

# ------------------------
# MYEGA temparature dependence of viscosity

# function obj_kin_myega(p, _)
#     # if p[1] < 0 #|| p[2] < 100 #|| p[3] < 100 
#     #     return NaN
#     # end
#     mk = k_f_myega(p)
#     ks = mk.(Tstops)
#     logks = log.(ustrip.(u"s^-1", ks))
#     return sum(abs2, logks .- klogs128)
# end
# objf_myega = OptimizationFunction(obj_kin_myega, AutoForwardDiff())
# objf_myega(p0_myega, nothing)

# # opt_kin_myega = solve(OP(objf_myega, p0_myega, lb=fill(0.0, 3), ub=fill(Inf, 3)), alg)
# opt_kin_myega = solve(OP(objf_myega, p0_myega), alg)

# -------------------------------------
# sigmoid temparature dependence of viscosity
# function err_kin_sigm(p,_)
#     ks = k_f_sigm(p).(Tstops)
#     err = log.(ustrip.(u"s^-1", ks)) - klogs128
# end

# objf_sigm = OptimizationFunction(obj_kin_sigm, AutoForwardDiff())
# obj_kin_sigm(p0_sigm,nothing)
# opt_sigm = solve(NLSP(err_kin_sigm, p0_sigm), LevenbergMarquardt(), show_trace=Val(true))

# # -------------------- 
# # Assume a beta distribution and see results
# gr()
# betaplot(sol1_sigm128, legend=:topright, legend_columns=1)
# # betaplot(sol8_sigm128, ylim=(1e-10, 1e35), xscale=:log10, yscale=:log10, legend=:bottomleft, legend_columns=1)
# plot!(xlabel="Pore volume [mm³]", ylabel="Probability Density")
# savefig(plotsdir("changingdist_linear.png"))

# -----------------------

# function namefitparam_arb(sol, label)
#     po = sol.prob.p
#     nt = (label = label,
#         N = po.N,
#         S = po.S,
#         C = po.C,
#         Rs = po.Rs,
#         k_f = po.k_f,
#     )
# end
# fits_sol = [sol1_MC1, sol8_MC8, sol1_MC18, sol1_MC128, sol1_myega128, sol1_sigm128]
# fits_lab = ["1", "8", "1&8", "1&2&8", "MYEGA 1&2&8", "Tanh 1&2&8"]
# # fits_sigm = Optim.minimizer.([opt12sigm, opt18sigm, opt128sigm]) 
# rep = Table(map(namefitparam_arb, fits_sol, fits_lab))

# # safesave(datadir("exp_pro", "fit12.html"), pretty_table(HTML, report, header = ["Fit", "N", "S", "k0", "E", "C", "Rs"]))
# safesave(datadir("exp_pro", "fit_arb.html"), pretty_table(HTML, rep, header = ["Fit", "N", "S", "C", "Rs", "kf"]))
# # pretty_table(report, header = ["Fit", "N", "S", "k0", "E", "C", "Rs"])

# ---------------------
# MCMC sampling
using Turing


@model function sigma_likelihood(pf)

    Nfac ~ Uniform(-6, 6)
    Sfac ~ Uniform(-6, 6)
    Cfac ~ Uniform(-6, 6)
    Rfac ~ Uniform(-6, 6)
    k0fac ~ Uniform(-6, 6)
    kffac ~ Uniform(-6, 6)
    K1fac ~ Uniform(-6, 6)
    K2fac ~ Uniform(-6, 6)

    # k_args = exp10.([kfac, Kfac, θfac]) .* fit_guess[5:7]
    allvals = [Nfac, Sfac, Cfac, Rfac, k0fac, kffac, K1fac, K2fac]

    if pf[1] isa PrimaryDryFit
        sqerr = obj_pd(allvals, (trans_sigma, pf...), tweight=1e-2)
    else
        sqerr = objn_pd(allvals, (trans_sigma, pf...), tweight=1e-2)
    end

    if !isnan(sqerr)
        Turing.@addlogprob! -sqerr
    else
    #     # @info "check" newpo.N.val newpo.S.val newpo.C.val newpo.Rs.val sol.t[end] sol[1,end]
        Turing.@addlogprob! -Inf
    end
end

case1_params = sigma_likelihood(make_mc_pass([1]))
case128_params = sigma_likelihood(make_mc_pass([1, 2, 8]))
    
smallalg = NUTS(10, 0.65)
init_single = [rand(8).-0.5]
smallsamp = @time sample(case1_params, smallalg, MCMCSerial(), 20, 1; initial_params=init_single)
nthreads = 4
init_set = [(rand(8).-0.5) for i in 1:nthreads]
smallsamp128 = @time sample(case128_params, smallalg, MCMCSerial(), 20, nthreads, initial_params = init_set)
# smallsamp128 = @time sample(case128_params, smallalg, MCMCThreads(), 20, nthreads, initial_params = init_set)
# @profview sample(case1_params, smallalg, MCMCSerial(), 20, 1)
# plot(smallsamp)

alg = NUTS(75, 0.65)
samp = @time sample(case1_params, alg, MCMCThreads(), 500, 4)
safesave(datadir("exp_pro", "case1_myega_sampling.jld2"), Dict("samp" => samp, "smallsamp" => smallsamp))
samp128 = @time sample(case128_params, alg, MCMCThreads(), 500, 4)
safesave(datadir("exp_pro", "case128_myega_sampling.jld2"), @strdict smallsamp128 samp128 )

# # safesave(datadir("exp_pro", "case1_myega_sampling.jld2"), Dict("samp" => samp, "smallsamp" => smallsamp))

# samp = load(datadir("exp_pro", "case1_myega_sampling_#1.jld2"))["samp"]

# plot(samp)

