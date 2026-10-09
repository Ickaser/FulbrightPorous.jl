fp_plot_defaults()

casenames = ["A3", "A5", "A6", "B"]
casecolors = [
    RGB(157/255, 14/255, 20/255),
    RGB(4/255, 73/255, 148/255),
    RGB(56/255, 160/255, 86/255),
    RGB(0, 0, 0),
]
readcaseT(name) = CSV.read(datadir("Kosir_et_al_2024", "$(name)_tT.csv"), Table, header=["t", "T"])
readcaseD(name) = CSV.read(datadir("Kosir_et_al_2024", "$(name)_Df.csv"), Table, header=["D", "f"])

data = Table(map(casenames, casecolors) do case, color
    tT = map(readcaseT.(case)) do ti
        (t = ti.t*u"hr", T = ti.T*u"°C")
    end
    Df = map(readcaseD.(case)) do di
        (D = di.D*u"μm", f = di.f*u"μm^-1")
    end
    (;name=case, tT, Df, color)
end)

# Re-plot pore distributions for presentation
begin
plot(u"μm", u"μm^-1", xlabel=L"pore size, $D_p$ [$\mu$m]", ylabel=L"prob. density, $f_D$ [$\mu$m$^{-1}$]", guidefontfamily="Computer Modern", unitformat=(x,y)->x)
for case in data
    plot!(case.Df.D, case.Df.f, label="Exp. $(case.name)", color=case.color, lw=4)
end
plot!(size=(300, 250), ylims=(0, 0.03), yticks=0:0.01:0.03,)
end
savefig(plotsdir("dextran_Df.svg"))

# Hard-coded for A3, A5, A6, B
t_end = [(18.0u"hr", 24.0u"hr"),
        (19.0u"hr", 26.0u"hr"),
        (20.0u"hr", 26.0u"hr"),
        (17.0u"hr", 26.0u"hr")]
csolid = [0.05u"g/mL", 0.15u"g/mL", 0.2u"g/mL", 0.05u"g/mL"]

begin
Ap, Av = π.*get_vial_radii("10R").^2 .|> u"cm^2"
Vfill = 3u"mL"
hf0 = Vfill/Ap
Tsh = RampedVariable([-40.0, -5]u"°C" .|> u"K", 1u"K/minute")
pch = RampedVariable(10.0u"Pa")
Kshf_g = ConstPhysProp(10.0u"W/m^2/K")
Rp_g = RpFormFit(0.1u"cm^2*Torr*hr/g", 30.0u"cm*Torr*hr/g", 0.1u"cm^-1")
ρsolution = 1.0u"g/mL"
cbase = 0.05u"g/mL"

d_typical = 50u"μm"
V_typical = 4/3*π*d_typical^3
N_typical = 1/V_typical
# Dummy values
N = ConstPhysProp(N_typical)
S = ConstPhysProp(5*V_typical) # = 5*V_typ^2*N_typ
C = 2e-1u"m^0.5"
Rs = 0.5u"cm^2*Torr*hr/g"
k0 = 0.0u"1/s"

po_base = ParamObjPikal(Rp_g, hf0, cbase, ρsolution, Kshf_g, Av, Ap, pch, Tsh) 
pomc_base = ParamObjMC([
    (hf0, cbase, ρsolution),
    (Kshf_g, Av, Ap),
    (pch, Tsh),
    (ConstPhysProp(k0), C, N, S, Rs)
])
end

cases = Table(data, (;t_end, csolid))

pdends = [13u"hr", 10.5u"hr", 14u"hr", 13.8u"hr"]
cases_pdfit = Table(map(cases, pdends) do case, pdend
    (; tT, t_end) = case
    pdfit = @df filter(x->x.t<pdend, tT) ExpFitData(:t .- :t[1], TfData(:T), EndTimeData(t_end .- :t[1]))
    po = @set po_base.csolid = case.csolid
    V = π/6*case.Df.D.^3
    μ1 = sum(case.Df.f .* V)/sum(case.Df.f) |> u"mL"
    μ0 =  Vfill / μ1
    N_sem = μ0 / Vfill
    # μ1 / Vfill
    μ2 = sum(case.Df.f .* V.^2)/sum(case.Df.f) |> u"mL^2"
    S_sem = μ2 / μ1
    @info "check" case.name μ2 μ1 cbrt(sqrt(μ2-μ1^2))|>u"μm" S_sem-1/N_sem
    pomc = @set pomc_base.csolid = case.csolid
    @reset pomc.N = ConstPhysProp(N_sem)
    @reset pomc.S = ConstPhysProp(S_sem)
    merge(case, (;pdfit, po, pomc, N_sem, S_sem))
end)

begin
pls = map(cases_pdfit) do case
    (; tT, name, t_end, pdfit) = case
    blankplot_hrC()
    # plot!(tT.t, tT.T, label="T")
    plot!(case.pdfit, label=["T for fit" L"t_\mathrm{end}"])
end
plot(pls...)
end

sol = solve(ODEProblem(po_base), Rodas3())
blankplot_hrC()
modconvtplot!(sol)
plot!(cases_pdfit.pdfit[1], )

# Fit to get Kshf
trans_KRp = KRp_transform_basic(Kshf_g.val, Rp_g.R0, Rp_g.A1, Rp_g.A2)
trans_Rp = Rp_transform_basic(Rp_g.R0, Rp_g.A1, Rp_g.A2)
shared_trans = as((separate = as(Vector, trans_Rp, 4),
      shared = K_transform_basic(Kshf_g.val),
      ))

objf_pd = OptimizationFunction((x,y)->obj_pd(x,y,tweight=1e-2), AutoForwardDiff())
objnf_pd = OptimizationFunction((x,y)->LyoPronto.objn_pd(x,y,tweight=1e-1), AutoForwardDiff())
@time objf_pd(fill(0.0, 4), (trans_KRp, po_base, cases_pdfit.pdfit[1]))
@time objnf_pd(fill(0.0, 13), (shared_trans, cases_pdfit.po, cases_pdfit.pdfit))

# Fit Kv to each separately
# KRp_opts = map(cases_pdfit) do case
#     po = @set po_base.csolid = case.csolid
#     @time opt_KRp = solve(OP(objf_pd, fill(0.0, 4), (trans_KRp, po, case.pdfit)), optalg)
# end
# Fit Kv to all at once
p0 = fill(0.1, 13); p0[1] = 2.8
p0[3] = log(1.5)
p0[9] = p0[12] = log(0.5)
@time opt_Kv = solve(OP(objnf_pd, p0, (shared_trans, cases_pdfit.po, cases_pdfit.pdfit)), optalg)

fitKv = transform(shared_trans, opt_Kv.u)

cases_Rp = Table(map(cases_pdfit, fitKv.separate) do case, fit
    po = setproperties(case.po, merge(fit, fitKv.shared))
    pomc = setproperties(case.pomc, fitKv.shared)
    hRp = calc_hRp_T(po, case.pdfit)
    return merge(case, (;po, pomc, hRp))
end)
psols = [solve(ODEProblem(po), Rodas3()) for po in cases_Rp.po]
# psols = LyoPronto.gen_nsol_pd(opt_Kv.u, shared_trans, cases_pdfit.po, cases_pdfit.pdfit);

begin
KRp_plots = map(cases_Rp, psols) do case, sol
    blankplot_hrC(legend=:bottomright, title=case.name)
    modconvtplot!(sol, c=:crimson)
    plot!(case.pdfit, c=case.color, nmarks=15)
    plot!(xlim=(-0.5, 21), ylim=(-45, -15))
end
plot(KRp_plots...)
end

# Ignore microcollapse kinetics in fit: treat as known at 0
porous_nt = FulbrightPores.porous_transform_bounded(Vfill)

NTV, STV, CTV, RsTV = map(key->getproperty(porous_nt, key), (:N, :S, :C, :Rs))
trans_cr = as((C=CTV, Rs=RsTV))
# Other transforms considered other fits:
# trans_nm = as((separate=as(Vector, as((N=NTV, S=STV)), 4),
#                 shared=as((C=CTV, Rs=RsTV, k_f=myega_nt.k_f)),))
# trans_ns_sem = as((separate=as(Vector, as((C=CTV,)), 4),
#                 shared=as((Rs=RsTV, k_f=sigma_nt.k_f)),))
# trans_ns = as((separate=as(Vector, as((N=NTV, S=STV)), 4),
#                 shared=as((C=CTV, Rs=RsTV, k_f=sigma_nt.k_f)),))

p0_cr = zeros(TransformVariables.dimension(trans_cr)) 
p0_cr = [0.5, 1.0]
blankplot_hrC()
modconvtplot!(LyoPronto.gen_nsol_pd(p0_cr, trans_cr, cases_Rp.pomc, cases_Rp.pdfit)..., c=[1 2 3 4])

@time opt_1234 = solve(OP(objnf_pd, p0_cr, (trans_cr, cases_Rp.pomc, cases_Rp.pdfit)), optalg)
@time opt_12 = solve(OP(objnf_pd, p0_cr, (trans_cr, cases_Rp.pomc[1:2], cases_Rp.pdfit[1:2])), optalg)

fit_cr = transform(trans_cr, opt_12.u)
# fit_cr = transform(trans_cr, opt_1234.u)

# fit_mc_sem = transform(trans_ns_sem, opt1234myegasem.u)
# fit_mc_sem_s = transform(trans_ns_sem_s, opt1234myegasems.u)
# all_pomc = [setproperties(pomc, merge(sep, fit_mc.shared)) for (pomc, sep) in zip(cases_Rp.pomc, fit_mc.separate)]
# all_pomc = [setproperties(pomc, merge(sep, fit_mc_sem.shared)) for (pomc, sep) in zip(cases_Rp.pomc, fit_mc_sem.separate)];
all_pomc = [setproperties(pomc, fit_cr) for pomc in cases_Rp.pomc];

sols_m = [solve(ODEProblem(pomc), Rodas3()) for pomc in all_pomc];
cases_both = Table(cases_Rp, (pomc=all_pomc,sol_mc=sols_m))

errs = Table(map(cases_both) do case
    PBM = case.sol_mc
    PP = solve(ODEProblem(case.po), Rodas3())
    errf = (model, tweight)->obj_expT(model, case.pdfit; tweight)
    PBM_T = sqrt(errf(PBM, 0))*u"K"
    PP_T = sqrt(errf(PP, 0))*u"K"
    t_end_est = ustrip(u"hr", sum(case.pdfit.t_end)/2)
    PBM_t = sqrt(errf(PBM, 1) - errf(PBM, 0))/t_end_est
    PP_t = sqrt(errf(PP, 1) - errf(PP, 0))/t_end_est
    return (;name=case.name, PBM_T, PP_T, PBM_t, PP_t)
end)

l = range(0u"cm", hf0, length=100)

begin
modcs = [2 4]
myega_plots = map(cases_both, psols) do case, solp
    blankplot_hrC(legend=:bottomright, )
    plot!(title=case.name)
    # plot!(Tsh, label=L"T_\mathrm{sh}", c=:black)
    plot!(case.pdfit, nmarks=10, c=case.color)

    modconvtplot!(solp, c=modcs[2], alpha=0.8, label="Pikal fit")
    mclab = (case.name ∈ ["A3", "A5"] ? "PBM fit" : "PBM eval")
    modconvtplot!(case.sol_mc, label=mclab, c=modcs[1], alpha=0.8)
    if case.name ∈ ["A5", "A6", "B"]
        plot!(ylabel="")
        plot!(yticks=(-40:5:-15, ""), left_margin=-13Plots.px)
    end
    plot!(xlim=(-0.5, 20.1), ylim=(-40, -14))
end

plot()
for (i, c) in enumerate(cases_both.color)
    vspan!([i-1,i], c=c, alpha=0.3, label="")
end
e1 = groupedbar!(errs.name, [errs.PP_T errs.PBM_T], color=modcs[[2 1]], label=["Pikal fit" "PBM fit"], ylabel=L"$T_\mathrm{f}$ error [$\degree$C]",unitformat=(x,y)->x, lw=1)
plot!(e1, legend=:bottom)

lRp = blankplot_lRp()
for (i, case) in enumerate(cases_both)
    hdT, RpT = case.hRp
    plot!(hdT, RpT, alpha=0.7, c=case.color, label="", seriestype=:samplemarkers, offset=i+3,step=3, m=:circle, la=0)
    mcrpplot!(case.sol_mc, c=case.color, lw=3, label="")
end
plot!([Inf], [Inf], c=:gray, lw=3, label="PBM")
scatter!([Inf], [Inf], c=:gray, label=L"T_\mathrm{f}(t)")
plot!(xlim=(0, 0.4), ylim=(0, 15), legend=:bottomright)
plot!(xlabel=L"dry layer thickness, $h_\mathrm{d}$ [cm]", ylabel=L"$R_p$ [cm Torr hr g$^{-1}$]", guidefontfamily="Computer Modern", unitformat=(x,y)->x)

plot(myega_plots..., e1, lRp, layout=@layout([grid(1,4); grid(1,2)]), size=(800,500))
plot!(subplot=1, left_margin=15Plots.px)
plot!(subplot=5, left_margin=15Plots.px)
plot!(bottom_margin=15Plots.px)
end
savefig(plotsdir("dextran_T_semC.svg"))
savefig(plotsdir("dextran_T_semC.pdf"))

begin
modcs = [:teal :crimson]
myega_plots = map(cases_both, psols) do case, solp
    @info "check" obj_expT(case.sol_mc, case.pdfit, tweight=0) obj_expT(case.sol_mc, case.pdfit, tweight=1)
    blankplot_hrC(legend=:bottomright, )
    plot!(title=case.name)
    # plot!(Tsh, label=L"T_\mathrm{sh}", c=:black)
    plot!(case.pdfit, c=case.color, nmarks=15)
    # vspan!([Inf], [Inf], c=case.color, alpha=0.4, label=L"$t_\mathrm{end}$,  on-off")

    modconvtplot!(solp, c=modcs[2], alpha=0.8, label="Pikal fit")
    modconvtplot!(case.sol_mc, label="PBM fit", c=modcs[1], alpha=0.8)
    # if case.name ∈ ["A3", "A5"]
    #     plot!(xlabel="")
    # end
    if case.name ∈ ["A5", "A6", "B"]
        plot!(ylabel="")
        plot!(yticks=(-45:5:-15, ""), left_margin=-5Plots.px)
    end
    plot!(xlim=(-0.5, 20.1), ylim=(-45, -14))
end

plot()
for (i, c) in enumerate(cases_both.color)
    vspan!([i-1,i], c=c, alpha=0.3, label="")
end
e1 = groupedbar!(errs.name, [errs.PP_T errs.PBM_T], color=modcs[[2 1]], label="", ylabel=L"$T_\mathrm{f}$ error [$\degree$C]",unitformat=(x,y)->x, lw=1)
# e2 = groupedbar(errs.name, [errs.PBM_t errs.PP_t], ylabel=L"$t$ rel. error", c=[2 3], lw=1)
plot(myega_plots..., e1, layout=@layout([grid(1,4) e{0.2w}]), size=(800,200))
plot!(subplot=1, left_margin=15Plots.px)
plot!(bottom_margin=15Plots.px)
end
# plot!(subplot=2, ylabel="hh")
savefig(plotsdir("dextran_T_semC_poster.svg"))
savefig(plotsdir("dextran_T_semC_poster.pdf"))

begin
blankplot_lRp()
for (i, case) in enumerate(cases_both)
    hdT, RpT = case.hRp
    plot!(hdT, RpT, alpha=0.7, c=case.color, label="", seriestype=:samplemarkers, offset=i+3,step=3, m=:o, lw=0)
    mcrpplot!(case.sol_mc, c=case.color, lw=3, label="")
end
plot!([Inf], [Inf], c=:gray, lw=3, label="PBM")
# plot!([Inf], [Inf], ls=:dash, c=:gray, alpha=0.6, lw=4, label="Pikal")
scatter!([Inf], [Inf], c=:gray, label=L"$R_p$ from T")
plot!(xlim=(0, 0.4), ylim=(0, 15), legend=:topleft)
plot!(size=(400, 400), xlabel=L"dry layer thickness, $h_\mathrm{d}$ [cm]", ylabel=L"$R_p$ [cm Torr hr g$^{-1}$]", guidefontfamily="Computer Modern", unitformat=(x,y)->x)
end
savefig(plotsdir("dextran_Rp.svg"))
savefig(plotsdir("dextran_Rp.pdf"))

begin
blankplot_lRp()
for (i, case) in enumerate(cases_both)
    hdT, RpT = case.hRp
    # plot!(subplot=4, l, case.Rp.(l), ls=:dash, c=parse(Int, case.name), label="3PP $(case.name)")
    # mcrpplot!(case.sol_mc, c=case.color, label="PBM $(case.name)")
    # plot!(l, case.po.Rp.(l), ls=:dash, c=case.color, label="3PP $(case.name)")
    # scatter!(hdT, RpT, alpha=0.9, c=case.color, label="Rp from T, $(case.name)", seriestype=:samplemarkers)
    plot!(hdT, RpT, alpha=0.7, c=case.color, label="", seriestype=:samplemarkers, offset=i+3,step=3, m=:o, lw=0)
    mcrpplot!(case.sol_mc, c=case.color, lw=3, label="")
    # plot!(l, case.po.Rp.(l), ls=:dash, lw=4, c=case.color, label="")
end
plot!([Inf], [Inf], c=:gray, lw=3, label="PBM")
# plot!([Inf], [Inf], ls=:dash, c=:gray, alpha=0.6, lw=4, label="Pikal")
scatter!([Inf], [Inf], c=:gray, label=L"$R_p$ from T")
plot!(xlim=(0, 0.4), ylim=(0, 15), legend=:topleft)
plot!(size=(250, 250), xlabel=L"dry layer thickness, $h_\mathrm{d}$ [cm]", ylabel=L"$R_p$ [cm Torr hr g$^{-1}$]", guidefontfamily="Computer Modern", unitformat=(x,y)->x)
end
savefig(plotsdir("dextran_Rp_poster.svg"))



# ns_SEM = Table(map(cases_Rp, fit_mc.separate, fit_mc_sem.separate) do case, fit, fitsem
#     # μ1 = π/6*sum(case.Df.f .* case.Df.D.^3)/sum(case.Df.f)
#     # N_sem = 1/μ1 |> u"mL^-1"
#     # S_sem = π^2/36 * sum(case.Df.f .* case.Df.D.^6)/sum(case.Df.f) / Vfill |> u"mL"
#     N_sem = case.pomc.N(0) |> u"μm^-3"
#     S_sem = case.pomc.S(0) |> u"μm^3"
#     # C_sem = fitsem.C

#     # N_mod = case.sol_mc[3,end] / Vfill
#     # S_mod = case.sol_mc[4,end]*u"mm"^6 / Vfill |> u"mL"
#     # N_fit = fit.N(0) |> u"mL^-1"
#     # S_fit = fit.S(0) |> u"mL"
#     # C_fit = fit_mc.shared.C
#     # C_shared = fit_mc_sem_s.C
#     # @info "check" N S sum(case.Df.f)
#     # return (;N_sem,N_mod, N_fit, S_sem, S_mod, S_fit, C_sem, C_fit)
#     # return (;case.name, N_sem, N_fit, S_sem, S_fit, C_fit, C_sem, C_shared)
# end)

begin
relevant_params = [:N_sem, :S_sem]
tab_filtered = hcat([getproperty(cases_pdfit, p) for p in relevant_params]...)
row_labels = LaTeXString.("{".*cases_pdfit.name .*"}")
column_labels = LaTeXString.(["{\$$(string(key)[1])\$, "*latexify(unit(getproperty(cases_pdfit[1], key)), fmt=SiunitxNumberFormatter())*"}" for key in relevant_params])
clipboard(latextabular(ustrip.(tab_filtered),
    side=row_labels, 
    adjustment=[Symbol("@{}"), :r, :S, :S, Symbol("@{}")],
    head=column_labels,
    fmt=SiunitxNumberFormatter(),
    latex=false,
    booktabs=true))
end