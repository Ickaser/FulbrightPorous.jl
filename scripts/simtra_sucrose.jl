fp_plot_defaults()
Latexify.set_default(labelformat=:square)

begin
# ------------------
# Load Kv function
# Computed in scripts/kv.jl, because it is used elsewhere too
Kshf = load(datadir("exp_pro", "simtra_kv.jld2"))["K_shf"]

# --------------------------
# Cycle 1 with original filenames, but labeled as case 4 in Table 8: use that numbering
dat4raw = CSV.read(datadir("exp_raw", "Simtra_Rp", "Suc_Cycle4.csv"), Table)
dat4 = map(simtra_rename, dat4raw)

t_end4 = identify_pd_end(dat4.t, dat4.pch_pir, Val(:onoff), window_width=391,  tmin=10u"hr",  tmax=85u"hr")

pch4 = RampedVariable(100u"mTorr")
# @df dat4 plot(:t, :pch_cm, ylim=(0, 200))
# @df dat4 plot(:t, :pch_pir, ylim=(0, 200))
# plot!(pch4)
# tendplot!(t_end4)

# There's a slight delay before shelf temperature ramps
Tsh4 = RampedVariable(uconvert.(u"K", [-40, -40, -25]u"°C"), [1u"K/minute", 0.5u"K/minute"], [23u"minute"])
# @df dat4 plot(:t, :Tsh, xlim=(0, 2))
# plot!(Tsh4)

fitdat4 = PrimaryDryFit(dat4.t, (
    filter(x->(!isnan(x) && x<-32.5u"°C"), dat4.T1),
    filter(!isnan, dat4.T2[dat4.t .< 60u"hr"]),
    # collect(filter(!isnan, dat4.T3)),
    ),
    t_end=t_end4
)
# plot(fitdat4)

# ------------------------------
# Cycle 6 according to publication naming

dat6raw = CSV.read(datadir("exp_raw", "Simtra_Rp", "Suc_Cycle6.csv"), Table)
dat6 = map(simtra_rename, dat6raw)
filter!(x->x.t .> 0u"hr", dat6)

# @df dat6 exptfplot!(:t, :T1, :T2, :T3, c=:red)
# plot!(ylim=(-40, -30))

t_end6 = identify_pd_end(dat6.t, dat6.pch_pir, Val(:onoff), window_width=391,  tmin=10u"hr",  tmax=120u"hr")
pch6 = RampedVariable(100u"mTorr")
# @df dat6 plot(:t, :pch_cm, ylim=(0, 200))
# @df dat6 plot(:t, :pch_pir, ylim=(0, 200))
# plot!(pch6)
# tendplot!(t_end6)

fitdat6 = @df dat6 PrimaryDryFit(:t, (
    # collect(filter(!isnan, :T1)),
    filter(!isnan, :T2[:t .< 85u"hr"]),
    filter(!isnan, :T3),),
    t_end=t_end6)
# plot(fitdat6)

Tsh6 = RampedVariable(uconvert.(u"K", [-45, -45, -30]u"°C"), [1u"K/minute", 0.5u"K/minute"], [23u"minute"])
# @df dat6 plot(:t, :Tsh, xlim=(0, 2))
# plot!(Tsh6)

# -------------------------------  
# Cycle 5
dat5raw = CSV.read(datadir("exp_raw", "Simtra_Rp", "Suc_Cycle5.csv"), Table)
dat5 = map(simtra_rename, dat5raw)

pch5 = RampedVariable(60u"mTorr")
t_end5 = identify_pd_end(dat5.t, dat5.pch_pir, Val(:onoff), window_width=391,  tmin=10u"hr",  tmax=68u"hr")
# @df dat5 plot(:t, :pch_cm, ylim=(0, 200))
# @df dat5 plot!(:t, :pch_pir, ylim=(0, 200))
# plot!(pch5)
# tendplot!(t_end5)

Tsh5 = RampedVariable(uconvert.(u"K", [-40, -25]u"°C"), 0.5u"K/minute")
# @df dat5 plot(:t, :Tsh, xlim=(0, 2))
# plot!(Tsh5)
fitdat5 = @df dat5 PrimaryDryFit(:t, (filter(!isnan, :T1), filter(!isnan, :T2), :T3[:t.<45u"hr"]), t_end=t_end5)
# plot(fitdat5)

# -----------------
# Cycle 5
dat2raw = CSV.read(datadir("exp_raw", "Simtra_Rp", "Suc_Cycle2.csv"), Table)
dat2 = map(simtra_rename, dat2raw)

pch2 = RampedVariable(50u"mTorr")
t_end2 = identify_pd_end(dat2.t, dat2.pch_pir, Val(:onoff), window_width=391,  tmin=10u"hr",  tmax=50u"hr")
# @df dat2 plot(:t, :pch_cm, ylim=(0, 200))
# @df dat2 plot!(:t, :pch_pir, ylim=(0, 200))
# plot!(pch2)
# tendplot!(t_end2)

Tsh2 = RampedVariable(uconvert.(u"K", [-40, -40, -20]u"°C"), [1u"K/minute", 0.5u"K/minute"], [12u"minute"])
# @df dat2 plot(:t, :Tsh_i, xlim=(0, 2))
# plot!(Tsh2)

fitdat2 = @df dat2 PrimaryDryFit(:t, (filter(!isnan, :T1), filter(!isnan, :T3), :T2[:t.<40u"hr"]), t_end=t_end2)
# plot(fitdat2)

# -----------------
# Cycle 1, high temperature
dat1raw = CSV.read(datadir("exp_raw", "Simtra_Rp", "Suc_Cycle1.csv"), Table)
dat1 = map(simtra_rename, dat1raw)
# @df dat1 exptfplot(:t, :T1, :T2, :T3)
# @df dat1 plot!(:t, :Tsh_i, xlim=(0, 20), ylim=(-40, -20))

pch1 = RampedVariable(100u"mTorr")
t_end1 = identify_pd_end(dat1.t, dat1.pch_pir, Val(:onoff), window_width=391,  tmin=10u"hr",  tmax=41u"hr")
# @df dat1 plot(:t, :pch_cm, ylim=(0, 200))
# @df dat1 plot(:t, :pch_pir, ylim=(0, 200))
# plot!(pch1)
# tendplot!(t_end1)

Tsh1 = RampedVariable(uconvert.(u"K", [-40, -40, -10]u"°C"), [1u"K/minute", 0.5u"K/minute"], [4u"minute"])
# @df dat1 plot(:t, :Tsh_i, xlim=(0, 2))
# plot!(Tsh1)

fitdat1 = @df dat1 PrimaryDryFit(:t, (
    filter(!isnan, :T1),
    filter(!isnan, :T2[:t .< 17u"hr"]),
    filter(!isnan, :T3[:t .< 17u"hr"]),),
    t_end=t_end1);
# plot(fitdat1)
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
pop4 = new_Tp(pop1, Tsh4, pch4)
pop5 = new_Tp(pop1, Tsh5, pch5)
pop6 = new_Tp(pop1, Tsh6, pch6)
pomc2 = new_Tp(pomc1, Tsh2, pch2)
pomc4 = new_Tp(pomc1, Tsh4, pch4)
pomc5 = new_Tp(pomc1, Tsh5, pch5)
pomc6 = new_Tp(pomc1, Tsh6, pch6)
dat1 = Table(dat1, (ld=FulbrightPores.calc_l(dat1, pop1.hf0)))
dat2 = Table(dat2, (ld=FulbrightPores.calc_l(dat2, pop2.hf0)))
dat4 = Table(dat4, (ld=FulbrightPores.calc_l(dat4, pop4.hf0)))
dat5 = Table(dat5, (ld=FulbrightPores.calc_l(dat5, pop5.hf0)))
dat6 = Table(dat6, (ld=FulbrightPores.calc_l(dat6, pop6.hf0)))
end

cases = Table([ (name="1", pop=pop1, pomc=pomc1, fitdat=fitdat1, dat=dat1),
                (name="2", pop=pop2, pomc=pomc2, fitdat=fitdat2, dat=dat2),
                (name="4", pop=pop4, pomc=pomc4, fitdat=fitdat4, dat=dat4),
                (name="5", pop=pop5, pomc=pomc5, fitdat=fitdat5, dat=dat5),
                (name="6", pop=pop6, pomc=pomc6, fitdat=fitdat6, dat=dat6),
                ])

trans_Rp = Rp_transform_basic(0.8u"cm^2*Torr*hr/g", 14.0u"cm*Torr*hr/g", 1e0u"1/cm")
obj_Rpf = OptimizationFunction((x,y)->obj_pd(x,y,tweight=1e-2), AutoForwardDiff(chunksize=3))
# obj_Rpf = OptimizationFunction(obj_pd, AutoForwardDiff(chunksize=3))
@time obj_Rpf([0.0, 1.0, 0.0], (trans_Rp, pop1, fitdat1))

# obj_npf = OptimizationFunction((x,y)->objn_pd(x,y,tweight=1e-2), AutoForwardDiff(chunksize=3))

single_Rps = map(cases) do case
    (;name, pop, fitdat, dat) = case
    u0 = [0., 0., 0.]
    @time opt = solve(OP(obj_Rpf, u0, (trans_Rp, pop, fitdat)), optalg)#, g_tol=1e-6)
    SciMLBase.successful_retcode(opt.retcode) || @warn "Optimization did not converge for Rp case $name"
    Rp = transform(trans_Rp, opt.u).Rp
    sol = gen_sol_pd(opt.u, trans_Rp, pop, fitdat)
    @info "analyis" name opt Rp
    hRps = map(x->calc_hRp_T(case.pop, case.fitdat, i=x), 1:length(case.fitdat.Tfs))
    return (opt=opt.u, sol_p=sol, Rp, hRps)
end;

cases_Rps = Table(cases, single_Rps)

pls = map(cases_Rps) do case 
    plT = blankplot_hrC()
    plot!(case.fitdat, nmarks=15, title=case.name)
    modconvtplot!(case.sol_p)
    plot!(legend = :bottomright)
    # plot!(case.pop.Tsh, tmax=case.fitdat.t_end, c=:black, label=L"T_\mathrm{sh}")
    return plT
end
plot(pls..., layout=(3,2), size=(800, 1000))

l = range(0u"cm", hf0, length=100)
begin
pls = map(cases_Rps) do case
    plr = blankplot_lRp()
    plot!(case.dat.ld, case.dat.Rp, label="TDLAS")
    for i in 1:length(case.fitdat.Tfs)
        hdT, RpT = case.hRps[i]
        plot!(hdT, RpT, c=i+1, label="Rp from TP$i")
    end
    plot!(l, case.Rp.(l), c=5, label="3PP fit directly to T")
    plot!(ylim=(0, 10), legend = :bottomright)
    plot!(title="\$p_{ch}=\\mathrm{$(case.pop.pch(Inf*u"s")|>u"mTorr")}\$,\n \$T_{sh}=$(float(case.pop.Tsh(Inf*u"s")|>u"°C"))\$  ")
    plot!(top_margin=20Plots.px)
    plot!(bottom_margin=15Plots.px)
    case.name ∈ ["1", "5"] && plot!(left_margin=20Plots.px)
    plot!(legend=:none)
    return plr
end
explot = plot(framestyle=:none, legendfontsize=12, legend_position=:inside)
plot!([Inf], [Inf], label="TDLAS", c=1)
plot!([Inf], [Inf], label="Rp from TP1", c=2)
plot!([Inf], [Inf], label="Rp from TP2", c=3)
plot!([Inf], [Inf], label="Rp from TP3", c=4)
plot!([Inf], [Inf], label="3PP fit directly to T", c=5)
plot(pls..., explot, layout=(2,3), size=(900, 500))
end
savefig(plotsdir("simtra_suc_Rp.png"))
savefig(plotsdir("simtra_suc_Rp.svg"))
savefig(plotsdir("simtra_suc_Rp.pdf"))

begin
blankplot_lRp()
colors = [
    RGB(230/255, 108/255, 9/255),
    :gray, 
    RGB(157/255, 187/255, 84/255), 
    :purple, 
    RGB(3/255, 176/255, 82/255),
]
map(cases_Rps, colors) do case, c
    # @df filter(x->x.Rp>0u"m/s", case.dat[begin:40:end]) plot!(:ld, :Rp; c, linealpha=0.5, label="$(case.name), TDLAS")
    # plot!(l, case.Rp.(l); seriestype=:samplemarkers, marker=:o, step=20, ma=0.5, c, ls=:dash, la=0.9, label="$(case.name), \$T\$ fit")
    pbmtdlas!(nothing, case.dat, case.Rp; hf0, ls=[:solid :dash], alpha = [0.7 0.9], marker=[:circle :none], label="", c=c)
    # @df filter(x->x.Rp>0u"m/s", case.dat[begin:40:end]) plot!(:ld, :Rp; c, linealpha=0.7, label="")
    # plot!(l, case.Rp.(l); seriestype=:samplemarkers, marker=:o, step=20, ma=0.5, c, ls=:dash, la=0.9, label="")
end
plot!(ylim=(0, 10), xlim=(0u"m", hf0*0.99), size=(350, 300), )
# plot!(legend = (0.9, 0.7), legend_columns=1, right_margin=30Plots.px)
plot!(legend = :bottomright, legend_columns=1)
plot!(ylabel=L"$R_p$")
scatter!([Inf], [Inf], label="TDLAS Estimation", c=:gray, marker=:circle, )
plot!([Inf], [Inf], c=:gray, label="Pikal Fit to T", ls=:dash)
end
savefig(plotsdir("simtra_suc_Rp_overlay.svg"))


# function make_Rp_pass(casenums)
#     subdat = filter(x->parse(Int, x.name)∈casenums, cases_Rps)
#     (subdat.pop, subdat.fitdat)
# end
function make_mc_pass(casenums)
    subdat = filter(x->parse(Int, x.name)∈casenums, cases)
    (subdat.pomc, subdat.fitdat)
end

# rpfit_cases = Table(cases = [(1,), (2,), (4,), (5,), (6,), (1,6)],
#         p0 = fill([0.5, -0.3, 0], 7),
#         trans = fill(trans_Rp, 7))
# rpfits = map(rpfit_cases) do case
#     sol = @time solve(OP(obj_npf, case.p0, (case.trans, make_Rp_pass(case.cases)...)), optalg;)
#     # SciMLBase.strip_solution(sol) # make this serialize nicer
# end

# -----------------------------
# Fit k by interpolating at four values
# interp_nt = (kl1 = TV.Identity(),
#             kl2 = TV.Identity(),
#             kl3 = TV.Identity(),
#             kl4 = TV.Identity())

const TV = TransformVariables
T_interps = [-45.0, -34, -31, -27]u"°C" .|> u"K"
T_interps_trans = as((;T=Constant(T_interps)))
k_logs_trans = as((;klog=as((kl1 = TV.Identity(),
            kl2 = TVExp(),
            kl3 = TVExp(),
            kl4 = TVExp()))))
interp_nt = as((;k_f = as(InterpK, merge(k_logs_trans, T_interps_trans))))
porous_nt = FulbrightPores.porous_transform_bounded(fillvol)
myega_nt = as((;k_f=as(MYEGA_k, as((k0 = TVScale(1e11u"s^-1") ∘ TVExp(),
            K1 = TVScale(5e4u"K") ∘ TVShift(0.001) ∘ TVLogistic() ∘ TVShift(TV.logit(0.05)),
            K2 = TVScale(6e2u"K") ∘ TVShift(0.01) ∘ TVLogistic() ,)))))
# plot(Tran, transform(myega_nt, zeros(3)).k_f.(Tran), yscale=:log10, label="MYEGA, def")
# sigma_nt = as((;k_f=as((k0 = TVScale(1e-15u"s^-1") ∘ TVExp(),
#             kf = TVScale(1e-8u"s^-1") ∘ TVExp(),
#             K1 = TVScale(0.8u"K") ∘ TVExp(),
#             K2 = TVScale(1.0u"K") ∘ TVShift(216.0) ∘ TVScale(45.0) ∘ TVLogistic() ∘TVScale(0.1) ,)
# )))
sigma_nt = FulbrightPores.sigmoid_transform_bounded(1e-8u"s^-1", 0.8u"K")
trans_interp = merge(porous_nt, interp_nt)
trans_myega = merge(porous_nt, myega_nt)
trans_sigma = merge(porous_nt, sigma_nt)

obj_nmc = OptimizationFunction((x,y)->objn_pd(x,y; badprms=FulbrightPores.badNS, tweight=1e-2), AutoForwardDiff(chunksize=4))


interp_guess = [-1.0, 1, 0, 0, -1, 0, 1, -2]
# plot(u"°C", u"s^-1")
# kineticplot!(transform(trans_interp, interp_guess),  label="Interp guess")
myega_guess = [0.1, 0.1, 0, 0, -10, 0, 0]
# kineticplot(transform(trans_myega, myega_guess), label="MYEGA guess")
sigma_guess = [-1.0, 1, 0, 0, 0, 0, 0, 0]
# kineticplot(transform(trans_sigma, sigma_guess), label="MYEGA guess")

mcfit_cases = Table(cases = repeat([(1,), (6,), (1, 6), (1, 4, 6)], 3),
        p0 = vcat(fill(interp_guess, 4), fill(myega_guess, 4), fill(sigma_guess, 4)),
        trans = vcat(fill(trans_interp, 4), fill(trans_myega, 4), fill(trans_sigma, 4)))
# debug
# mcfits = map(mcfit_cases[1:1]) do case
#     po, fitdat = make_mc_pass(case.cases)
#     po = po[1]
#     fitdat = fitdat[1]
#     newp = setproperties(po, transform(case.trans, case.p0))
#     prob = ODEProblem(newp)
#     solve(prob, odealg_chunk4)
# end

mcfits = map(mcfit_cases) do case
    sol = @time solve(OP(obj_nmc, case.p0, (case.trans, make_mc_pass(case.cases)...)), optalg;)
    if !SciMLBase.successful_retcode(sol.retcode)
        @warn "Optimization did not converge for cases $(case.cases), transfrom $(case.trans)"
    end
    # @time solve(OP(obj_nmc, case.p0, (case.trans, make_mc_pass(case.cases)...)), optalg; show_trace=true)
    sol.u
    # SciMLBase.strip_solution(sol) # make this serialize nicer
end

begin
plot(u"°C", u"s^-1")
marks = Dict(3=>:circle, 7=>:square, 11=>:utriangle)
for i in [3,7,11]
    k = transform(mcfit_cases[i].trans, mcfits[i]).k_f
    kineticplot!(k, markershape=marks[i], label="$(nameof(typeof(k))) fit to $(mcfit_cases[i].cases)")
end
# kineticplot!(transform(trans_interp, interp_guess),  label="Interp guess")
plot!(ylim=(1e-20, 2e-5), legend=:bottomright)
end



function get_sol_mc(casenum, fitnum)
    case = filter(x->parse(Int, x.name) == casenum, cases_Rps)[1]
    sol = gen_sol_pd(mcfits[fitnum], mcfit_cases[fitnum].trans, case.pomc)
end
# function get_sol_Rp(casenum, fitnum)
#     case = filter(x->parse(Int, x.name) == casenum, cases_Rps)[1]
#     sol = gen_sol_pd(rpfits[fitnum], rpfit_cases[fitnum].trans, case.pop)
# end

# Choose a favored fit for plotting, etc.
# fitnum = 3 # 2&8, interp
# fitnum = 7 # 2&8, MYEGA
fitnum = 11 # 2&8, sigmoid

bestfit = transform(mcfit_cases[fitnum].trans, mcfits[fitnum])
@info "Chosen fit" bestfit
bestfit.N.val |> u"μm^-3"
bestfit.S.val/fillvol |> NoUnits
bestfit.C |> u"m^0.5"

errs = Table(map(cases_Rps[[1, 2, 3, 4, 5]]) do case
    i = parse(Int, case.name)
    PBM = get_sol_mc(i, fitnum)
    # PP = get_sol_Rp(i, 6)
    PPi = case.sol_p
    errf = (model, tweight)->obj_expT(model, case.fitdat; tweight)
    PBM_T = sqrt(errf(PBM, 0))*u"K"
    # PP_T = sqrt(errf(PP, 0))*u"K"
    PPi_T = sqrt(errf(PPi, 0))*u"K"
    t_end_est = ustrip(u"hr", sum(case.fitdat.t_end)/2)
    PBM_t = sqrt(errf(PBM, 1) - errf(PBM, 0))/t_end_est
    @info "check" PBM_t PBM.t[end] t_end_est
    # PP_t = sqrt(errf(PP, 1) - errf(PP, 0))/t_end_est
    PPi_t = sqrt(errf(PPi, 1) - errf(PPi, 0))/t_end_est
    pch = case.pop.pch(Inf*u"s") |> u"mTorr"
    Tsh = case.pop.Tsh(Inf*u"s") |> u"°C" |> float
    return (;name=case.name, PBM_T, PBM_t, PPi_T, PPi_t, pch, Tsh)
end)

# Save results to disk
safesave(datadir("exp_pro", "simtra_sucrose_fitdata.jld2"), Dict("fitdats"=>cases.fitdat))
safesave(datadir("exp_pro", "simtra_sucrose_fitresults.jld2"), @strdict errs mcfit_cases mcfits cases_Rps)

# Load results from disk for replotting without recomputing
fitdats = load(datadir("exp_pro", "simtra_sucrose_fitdata.jld2"))["fitdats"]
results = load(datadir("exp_pro", "simtra_sucrose_fitresults.jld2"))
@unpack errs, mcfit_cases, mcfits, cases_Rps = results


begin
modcs = [2, 3]
permute = [1, 5, 2, 3, 4]
labs = ["PBM fit to $(join(mcfit_cases[fitnum].cases, ", ", " & "))" ]
caselabs = ["p=$(ustrip(p))\nT=$(Int(ustrip(T)))" for (p, T) in zip(errs.pch, errs.Tsh)]
# Terr = groupedbar(errs.name, 
plot()
vspan!([0  1; 1  2], color=:gray, alpha=0.3, label="")
Terr = groupedbar!(caselabs,
    [errs.PPi_T[permute] errs.PBM_T[permute] ], c=modcs[[2 1]], 
    label=labs, lw=1)
plot!(xlabel="Case", ylabel=L"$T_\mathrm{f}$ error", yunitformat=:square, legend=:topleft)
plot!(xgrid=false, xminorticks=2, xtickdir=:none, xminortickdir=:none, xminorgrid=true, minorgridalpha=0.1, ygridalpha=0.2)
# plot!(xticks=(1:5, ""))
# terr = groupedbar(errs.name, 
plot()
vspan!([0  1; 1  2], color=:gray, alpha=0.3, label="")
terr = groupedbar!(caselabs,
    # [errs.PBM_t errs.PPi_t errs.PP_t ], c=[2 3 4], 
    [errs.PPi_t[permute] errs.PBM_t[permute] ], c=modcs[[2 1]], 
    label=labs, lw=1)
plot!(xlabel="Case", ylabel=L"$t$ relative error", yunitformat=:square, legend=:topleft)
plot!(xgrid=false, xminorticks=2, xtickdir=:none, xminortickdir=:none, xminorgrid=true, minorgridalpha=0.1, ygridalpha=0.2)
plot!(yticks=(0:0.1:0.2))

err_plot = plot(Terr, terr, layout=(2,1), size=(600,600))

pls = map(cases_Rps[permute], ) do case
    i = parse(Int, case.name)
    pl = blankplot_hrC(legend=:bottomright)
    # plot!(title=case.name)
    lab_pch = case.pop.pch(0)
    lab_Tsh = ustrip(u"degC",case.pop.Tsh(Inf*u"hr")) |> Int
    title =L"p_\mathrm{ch}=%$(latexify(lab_pch, env=:raw))"*"\n"*L"$T_\mathrm{sh}=%$(lab_Tsh) ^\circ C$"
    # lab_Tsh = round(typeof(1u"degC"), case.pop.Tsh(Inf*u"hr")) 
    # title = L"$p_\mathrm{ch}="*latexify(lab_pch, env=:raw)*L"$"*"\n"*L"$T_\mathrm{sh}="*latexify(lab_Tsh, env=:raw)*L"$"
    
    plot!(;top_margin=20Plots.px, title)
    nmarks = Int(ustrip(u"hr", case.fitdat.t_end[1]) ÷ 10)
    # scatter!([55], [-39], label="", ms=10, msw=0, mc = casec)
    plot!(case.fitdat, nmarks=nmarks)
    # modconvtplot!(get_sol_mc(i, fitnum), label="PBM, fit to $(mcfit_cases[fitnum].cases)")
    # modconvtplot!(case.sol_p, c=:green, alpha=0.6, label="3PP, fit to $(case.name)")
    # modconvtplot!(get_sol_Rp(i, 6), c=:purple, label="3PP, fit to $(rpfit_cases[6].cases)")
    # modconvtplot!(get_sol_Rp(i, 6), label="Pikal individual fit", c=modcs[2])
    modconvtplot!(case.sol_p, label="Pikal fit", c=modcs[2])
    modconvtplot!(get_sol_mc(i, fitnum), label="PBM fit to hi-lo", c=modcs[1])
    plot!(xlim=(-1, 120), ylim=(-40, -28), yticks=[-40, -35, -32, -29]) 
    plot!(legend=:none)
    i ∈ [2, 4, 6] && plot!(ylabel="")
    # i ∈ [1, 2, 3] && plot!(xlabel="")
    i ∈ [1, 6] && plot!(bg_legend = :gray95)
    if i ∈ [1] 
        plot!(left_margin=20Plots.px)
    else
        plot!(ylabel="", left_margin=-5Plots.px, yticks=[-40, -36, -32, -29])
    end
    i ∈ [5, 6] && plot!(bottom_margin=20Plots.px)
    return pl
end
plot!(pls[5], legend=(-0.3, -0.5), legendfontsize=12)
plot!(err_plot, legend=:none, bottom_margin=20Plots.px, left_margin=20Plots.px)
# plot!(err_plot[1], xlabel="")
# plot!(err_plot[2], bottom_margin=20Plots.px)
# plot(pls..., err_plot, layout=(2,3), size=(1000, 500))
plot(pls..., plot(err_plot[1]), plot(err_plot[2], left_margin=10Plots.px), 
    layout=@layout([grid(1,5){0.6h} ; b{0.3w} a{0.3w} _{0.4w} _{0.1w}]), 
    size=(1000, 440))
end
savefig(plotsdir("simtra_suc_fitval.svg"))
# no PDF because takes a little manual tweaking for a good legend
function condlabel(pop)
    lab_pch = pop.pch(Inf*u"s") |> u"mTorr"
    lab_Tsh = ustrip(u"degC",pop.Tsh(Inf*u"s")) |> Int
    return L"p_\mathrm{ch}=%$(latexify(lab_pch, env=:raw))"*"\n"*L"$T_\mathrm{sh}=%$(lab_Tsh) ^\circ C$"
end

begin
modcs = [:teal :crimson]
permute = [1, 5, 2, 3, 4]
labs = ["PBM fit to $(join(mcfit_cases[fitnum].cases, ", ", " & "))" ]
caselabs = ["p=$(ustrip(p))\nT=$(Int(ustrip(T)))" for (p, T) in zip(errs.pch, errs.Tsh)]
# Terr = groupedbar(errs.name, 
plot()
vspan!([0  1; 1  2], color=:gray, alpha=0.3, label="")
Terr = groupedbar!(caselabs,
    [errs.PPi_T[permute] errs.PBM_T[permute] ], c=modcs[[2 1]], 
    label=labs, lw=1)
plot!(xlabel="Case", ylabel=L"$T_\mathrm{f}$ error", yunitformat=:square, legend=:topleft)
plot!(xgrid=false, xminorticks=2, xtickdir=:none, xminortickdir=:none, xminorgrid=true, minorgridalpha=0.1, ygridalpha=0.2)
# plot!(xticks=(1:5, ""))
# terr = groupedbar(errs.name, 
plot()
vspan!([0  1; 1  2], color=:gray, alpha=0.3, label="")
terr = groupedbar!(caselabs,
    # [errs.PBM_t errs.PPi_t errs.PP_t ], c=[2 3 4], 
    [errs.PPi_t[permute] errs.PBM_t[permute] ], c=modcs[[2 1]], 
    label=labs, lw=1)
plot!(xlabel="Case", ylabel=L"$t$ relative error", yunitformat=:square, legend=:topleft)
plot!(xgrid=false, xminorticks=2, xtickdir=:none, xminortickdir=:none, xminorgrid=true, minorgridalpha=0.1, ygridalpha=0.2)
plot!(yticks=(0:0.1:0.2))

err_plot = plot(Terr, terr, layout=(2,1), size=(600,600))

pls = map(cases_Rps[permute], colors[permute]) do case, casec
    i = parse(Int, case.name)
    pl = blankplot_hrC(legend=:bottomright)
    # plot!(title=case.name)
    title = condlabel(case.pop)
    # lab_Tsh = round(typeof(1u"degC"), case.pop.Tsh(Inf*u"hr")) 
    # title = L"$p_\mathrm{ch}="*latexify(lab_pch, env=:raw)*L"$"*"\n"*L"$T_\mathrm{sh}="*latexify(lab_Tsh, env=:raw)*L"$"
    
    plot!(;top_margin=20Plots.px, title)
    nmarks = Int(ustrip(u"hr", case.fitdat.t_end[1]) ÷ 10)
    # scatter!([55], [-39], label="", ms=10, msw=0, mc = casec)
    plot!(case.fitdat, nmarks=nmarks, c=casec)
    # modconvtplot!(get_sol_mc(i, fitnum), label="PBM, fit to $(mcfit_cases[fitnum].cases)")
    # modconvtplot!(case.sol_p, c=:green, alpha=0.6, label="3PP, fit to $(case.name)")
    # modconvtplot!(get_sol_Rp(i, 6), c=:purple, label="3PP, fit to $(rpfit_cases[6].cases)")
    # modconvtplot!(get_sol_Rp(i, 6), label="Pikal individual fit", c=modcs[2])
    modconvtplot!(case.sol_p, label="Pikal fit", c=modcs[2])
    modconvtplot!(get_sol_mc(i, fitnum), label="PBM fit to hi-lo", c=modcs[1])
    plot!(xlim=(-1, 120), ylim=(-45, -28), yticks=[-45, -40, -35, -32, -29]) 
    i ∈ [2, 4, 6] && plot!(ylabel="")
    # i ∈ [1, 2, 3] && plot!(xlabel="")
    i ∈ [1, 6] && plot!(bg_legend = :gray95)
    if i ∈ [1] 
        plot!(left_margin=20Plots.px)
    else
        plot!(ylabel="", left_margin=-5Plots.px)
    end
    i ∈ [5, 6] && plot!(bottom_margin=20Plots.px)
    return pl
end
plot!(err_plot, legend=:none, bottom_margin=20Plots.px, left_margin=20Plots.px)
# plot!(err_plot[1], xlabel="")
# plot!(err_plot[2], bottom_margin=20Plots.px)
# plot(pls..., err_plot, layout=(2,3), size=(1000, 500))
plot(pls..., plot(err_plot[1]), plot(err_plot[2], left_margin=10Plots.px), 
    layout=@layout([grid(1,5){0.6h} ; b{0.3w} a{0.3w} _{0.4w} _{0.1w}]), 
    size=(1000, 440))
end
savefig(plotsdir("simtra_suc_fitval_poster.svg"))


begin
porecolors = colormap("oranges")[[40, 55, 70, 80, 90]] # Chosen from among 100 colors
markers = [:pentagon, :diamond, :utriangle, :rect, :circle]
plot!(palette=:default)
rpattr = (seriestype=:samplemarkers, linestyle=:dash)
labs = ["PBM fit to $(join(mcfit_cases[fitnum].cases, ", ", " & "))" "3PP individual fits " ]
plot(layout=@layout([grid(3,1) a{0.4w}]), unitformat=:square, size=(700,400))
for (i, casenum) in enumerate(reverse([1, 4, 2, 5, 6]))
    # poreevolplot!(get_sol_mc(casenum, fitnum), label="Case $casenum", c=porecolors[i], kind=:stats, marker=markers[i])
    sol = get_sol_mc(casenum, fitnum)
    poreevolplot!(sol, label="\n"*condlabel(sol.prob.p), c=porecolors[i], kind=:stats, marker=markers[i], lw=2)
end
# poreevolplot!(get_sol_mc(1, fitnum), label="PBM $(join(mcfit_cases[fitnum].cases, ", ", " & ")), Cycle 1", c=2, kind=:stats, marker=:circle)
# poreevolplot!(get_sol_mc(2, fitnum), label="PBM $(join(mcfit_cases[fitnum].cases, ", ", " & ")), Cycle 2", c=2, kind=:stats, marker=:square)
# poreevolplot!(get_sol_mc(6, fitnum), label="PBM $(join(mcfit_cases[fitnum].cases, ", ", " & ")), Cycle 6", c=2, kind=:stats, marker=:utriangle)
plot!(subplot=2, ylim=(0, Inf))
plot!(subplot=3, ylim=(0, Inf))
plot!(subplot=4, ylabel="R_p", yunitformat=latexify)
plot!(subplot=4, legend=(1.25, 0.7), right_margin=130Plots.px)
# plot!(l, cases_Rps[1].Rp.(l), subplot=4, label="3PP 1"; c=3, marker=:circle, rpattr...)
# plot!(l, cases_Rps[2].Rp.(l), subplot=4, label="3PP 2"; c=3, marker=:square, rpattr...)
# plot!(l, cases_Rps[5].Rp.(l), subplot=4, label="3PP 6"; c=3, marker=:utriangle, rpattr...)
end
savefig(plotsdir("simtra_suc_poreevol.svg"))
# Manually export to PDF after fixing markers in legend

begin
plot(u"°C", u"s^-1", xlabel="T\\;", ylabel="k \\, ", unitformat=latexify, guidefontsize=15)
k = transform(mcfit_cases[3].trans, mcfits[3]).k_f
# lab = "$(nameof(typeof(k))) fit to $(mcfit_cases[fitnum].cases)"
lab = "Fit 1 to cycles"
kineticplot!(k, markershape=:circle, label=lab)
k = transform(mcfit_cases[11].trans, mcfits[11]).k_f
# lab = "$(nameof(typeof(k))) fit to $(mcfit_cases[fitnum].cases)"
lab = "Fit 2 to cycles"
kineticplot!(k, markeralpha=0, label=lab)
Tran = -45u"°C":2u"K":-23u"°C"
k_arrhenius = T->1e-7u"s^-1" * exp(-2000u"K"/T) # Near-maximum slope for Arrhenius
plot!(Tran, k_arrhenius.(Tran), label="Arrhenius with max slope", c=:gray, marker=:none, ls=:dash)
# kineticplot!(transform(trans_interp, interp_guess),  label="Interp guess")
# plot!(ylim=(1e-20, 2e-0))
plot!(legend=:topleft, size=(500,250), bottom_margin=20Plots.px)
end
savefig(plotsdir("simtra_suc_kinetics.svg"))

# -----------------
# Inspect fit results in Rp space rather than temperature space

begin
locs = [(0.1, 6), (0.2, 2), (0.2, 2), (0.2, 2), (0.2, 2)]
Tps = [-29.5, -33.3, -32.5, -33.2, -33.6]u"°C"
permute = [1, 3, 2, 4, 5]
pls = map(cases_Rps[permute], locs, Tps[permute]) do case, loc, Tp
    i = parse(Int, case.name)
    pl = blankplot_lRp(ylims=(0, 11), xlims=(0, 1.3), legend=:topleft)
    plot!([Inf], [Inf], label=L"$T_\mathrm{f}(t)$", c=3)
    pbmtdlas!(get_sol_mc(i, fitnum), case.dat, nothing; label=["TDLAS" L"PBM fit to $T_\mathrm{f}(t)$"])
    for j in 1:length(case.fitdat.Tfs)
        scatter!(calc_hRp_T(case.pop, case.fitdat; i=j)..., label="",c=3, markeralpha=0.02, ms=2, msw=0)
    end
    plot!(title=condlabel(case.pop), top_margin=25Plots.px)
    # Tfmax = maximum([maximum(Tfi) for Tfi in case.fitdat.Tfs])
    # annotate!(loc..., text("\$T_\\mathrm{f,max}=$(float(Tfmax|>u"°C"))\$", :left, 12))
    annotate!(loc..., text("\$T_\\mathrm{f}=$(float(Tp|>u"°C"))\$", :left, 12))
    # i ∈ [1, 2] && plot!(xlabel="")#, xformatter=:none)
    i ∈ [1] && plot!(left_margin=20Plots.px, yticks=(0:2:12))
    i ∈ [2, 4, 5, 6] && plot!(ylabel="", left_margin=-10Plots.px, yticks=(0:2:12, []), legend=:none)
    pl
end
plot(pls..., layout=(1,5), size=(1000,300), bottom_margin=20Plots.px)
end
savefig(plotsdir("simtra_suc_Rp_val.svg"))
savefig(plotsdir("simtra_suc_Rp_val.pdf"))

begin
pl = blankplot_lRp(ylim=(0, 10), xlim=(0u"cm", hf0*0.99))
for (i, case, c) in zip([1,2,4,5,6], cases_Rps, colors)
    pbmtdlas!(get_sol_mc(i, fitnum), case.dat, nothing; c, label="", ls=[:solid :solid], marker = [:circle :none], alpha=[0.7 0.9])
end
scatter!([Inf], [Inf], label="TDLAS Estimation", c=:gray, marker=:circle, )
plot!([Inf], [Inf], c=:gray, label="PBM Fit to hi-lo T", ls=:solid)
plot!(size=(350, 300))
end
savefig(plotsdir("simtra_suc_Rp_PBM_overlay.svg"))

# ------------------------
# Knudsen numbers
FulbrightPores.calc_Kn(get_sol_mc(6, fitnum))
FulbrightPores.calc_Kn(get_sol_mc(1, fitnum))
Table(vcat(map([1,2,4,5,6]) do i
    sol = get_sol_mc(i, fitnum)
    Kn_tab = FulbrightPores.calc_Kn(sol)[[begin+1, end]]
end...))




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

# # ---------------------
# # MCMC sampling
# using Turing


# @model function sigma_likelihood(pf)

#     Nv ~ LogUniform(1e-6, 1e0) # μm^-3
#     Sv ~ LogUniform(1e-16, 1e0) # cm^3
#     Cv ~ LogUniform(1e-5, 1e0) # m^(1//2)
#     Rv ~ LogUniform(1e-10, 1e0) # cm^2*Torr*hr/g
#     k0v ~ LogUniform(1e-20, 1e-10) # s^-1
#     kfv ~ LogUniform(1e-20, 1e-5) # s^-1
#     K1v ~ Uniform(0.1, 10.0) # K
#     K2v ~ Uniform(216, 270) # K


#     # k_args = exp10.([kfac, Kfac, θfac]) .* fit_guess[5:7]
#     # allvals = [Nfac, Sfac, Cfac, Rfac, k0fac, kffac, K1fac, K2fac]
#     newps = (;
#         N = ConstPhysProp(Nv * u"μm^-3"),
#         S = ConstPhysProp(Sv * u"cm^3"),
#         C = Cv * u"m^0.5",
#         Rs = Rv * u"cm^2*Torr*hr/g",
#         k_f = Sigmoid_k(
#             k0v * u"s^-1",
#             kfv * u"s^-1",
#             K1v * u"K",
#             K2v * u"K",
#         )
#     )
    

#     # if pf[2] isa PrimaryDryFit
#     #     p = setproperties(pf[1], newps)
#     #     sol = solve(ODEProblem(p), odealg_chunk4)
#     #     sqerr = obj_expT(sol, pf[2], tweight=1e-2)
#     # else
#         sqerr = mapreduce(+, pf...) do pomc, fitdat
#             p = setproperties(pomc, newps)
#             @info "check1" p.N
#             sol = solve(ODEProblem(p), odealg_chunk4)
#             @info "check2" sol[end] sol.retcode
#             obj_expT(sol, fitdat, tweight=1e-2)
#         end
#     # end

#     if !isnan(sqerr)
#         Turing.@addlogprob! -sqerr
#     else
#     #     # @info "check" newpo.N.val newpo.S.val newpo.C.val newpo.Rs.val sol.t[end] sol[1,end]
#         Turing.@addlogprob! -Inf
#     end
#     return newps
# end

# case1_params = sigma_likelihood(make_mc_pass([1]))
# case128_params = sigma_likelihood(make_mc_pass([1, 2, 8]))

    
# smallalg = NUTS(10, 0.65)
# # guess_params = (Nv = 1e-4, Sv = 1e-8, Cv = 0.1, Rv = 1e-2, k0v = 1e-14, kfv = 1e-10, K1v = 1.0, K2v = 250.0)
# # init_single = [InitFromParams(guess_params)]
# # Turing.DynamicPPL.TestUtils.AD.run_ad(case1_params, AutoForwardDiff(); test=false, benchmark=false)
# smallsamp = @time sample(case1_params, smallalg, MCMCSerial(), 20, 1; initial_params=init_single, adtype=AutoForwardDiff())

# bad_params = (Nv = 0.00025055, Sv = 7.344991e-8, Cv = 0.02147, Rv = 3.96238e-10, k0v = 8.38e-16, kfv = 1.444e-9, K1v = 0.3905, K2v = 239.8)
# init_single = [InitFromParams(bad_params)]
# pomc_bad = setproperties(pomc1, (case1_params | bad_params)())
# sol_bad = solve(ODEProblem(pomc_bad), odealg_chunk4)





# nthreads = 4
# init_set = [(rand(8).-0.5) for i in 1:nthreads]
# smallsamp128 = @time sample(case128_params, smallalg, MCMCSerial(), 20, nthreads, initial_params = init_set)
# # smallsamp128 = @time sample(case128_params, smallalg, MCMCThreads(), 20, nthreads, initial_params = init_set)
# # @profview sample(case1_params, smallalg, MCMCSerial(), 20, 1)
# # plot(smallsamp)

# alg = NUTS(75, 0.65)
# samp = @time sample(case1_params, alg, MCMCThreads(), 500, 4)
# safesave(datadir("exp_pro", "case1_myega_sampling.jld2"), Dict("samp" => samp, "smallsamp" => smallsamp))
# samp128 = @time sample(case128_params, alg, MCMCThreads(), 500, 4)
# safesave(datadir("exp_pro", "case128_myega_sampling.jld2"), @strdict smallsamp128 samp128 )

# # # safesave(datadir("exp_pro", "case1_myega_sampling.jld2"), Dict("samp" => samp, "smallsamp" => smallsamp))

# # samp = load(datadir("exp_pro", "case1_myega_sampling_#1.jld2"))["samp"]

# using Turing
# samp = load(datadir("exp_pro", "case1_myega_sampling.jld2"))["samp"]

# plot(samp, xlim=(-6, 6))
# samp[:Nfac]
# samp[2]

# # plot(samp)

