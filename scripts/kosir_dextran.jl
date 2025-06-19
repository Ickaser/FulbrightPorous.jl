
const OP = OptimizationProblem
using LineSearches
optalg = BFGS(linesearch=LineSearches.BackTracking())
default(:linewidth, 3)
default(:markersize, 5)
using UnitfulLatexify

casenames = ["A3", "A5", "A6", "B"]
readcaseT(name) = CSV.read(datadir("Deck_et_al_2024", "$(name)_tT.csv"), Table, header=["t", "T"])
readcaseD(name) = CSV.read(datadir("Deck_et_al_2024", "$(name)_Df.csv"), Table, header=["D", "f"])


data = Table(map(casenames) do case
    tT = map(readcaseT.(case)) do ti
        (t = ti.t*u"hr", T = ti.T*u"°C")
    end
    Df = map(readcaseD.(case)) do di
        (D = di.D*u"μm", f = di.f*u"μm^-1")
    end
    (;name=case, tT, Df)
end)

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
Rp_g = RpFormFit(0.1u"cm^2*Torr*hr/g", 14.0u"cm*Torr*hr/g", 0.1u"cm^-1")
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
k0 = 1e-14u"1/s"

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
    pdfit = @df filter(x->x.t<pdend, tT) PrimaryDryFit(:t .- :t[1], :T, t_end .- :t[1])
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
    merge(case, (;pdfit, po, pomc))
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
@time opt_Kv = solve(OP(objnf_pd, fill(0.0, 13), (shared_trans, cases_pdfit.po, cases_pdfit.pdfit)), optalg)

fitKv = transform(shared_trans, opt_Kv.u)

# cases_Rp = Table(cases_pdfit)
# cases_Rp.po .= [setproperties(po, merge(Rp, fitKv.shared)) for (po, Rp) in zip(cases_pdfit.po, fitKv.separate)]
# cases_Rp.pomc .= [@set pomc.Kshf = fitKv.shared.Kshf for pomc in cases_pdfit.pomc]
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
    modconvtplot!(sol)
    plot!(case.pdfit)
    plot!(xlim=(-0.5, 21), ylim=(-45, -15))
end
plot(KRp_plots...)
end

# function make_mcpass(casenums)
#     subdat = filter(x->parse(Int, x.name)∈casenums, cases)
#     (subdat.pomc, subdat.pdfit)
# end

# pass123 = make_mcpass([1, 2, 3])
# pass45 = make_mcpass([4, 5])

NTV = ConstWrapTV() ∘ TVScale(4e7u"cm^-3") ∘ TVExp()
STV = ConstWrapTV() ∘ TVScale(2e-4u"cm^3") ∘ TVExp()
CTV = TVScale(1.2u"m^0.5") ∘ TVExp()
RsTV = TVScale(1.5u"cm^2*Torr*hr/g") ∘ TVShift(0.01) ∘ TVLogistic()
porous_nt = (N = NTV, S = STV, C = CTV, Rs = RsTV,) # If I could enforce S > ϕ/N, would be great...
myega_nt = (k0 = TVScale(1e42u"s^-1") ∘ TVExp(),
            K1 = TVScale(2.0e4u"K") ∘ TVExp(),
            K2 = TVScale(1e2u"K") ∘ TVExp() ,)
# plot(Tran, FulbrightPores.k_f_nt(transform(as(myega_nt), zeros(3))).(Tran), yscale=:log10, label="MYEGA, def")
sigma_nt = (k0 = TVScale(1e-15u"s^-1") ∘ TVExp(),
            kf = TVScale(1e-8u"s^-1") ∘ TVExp(),
            K1 = TVScale(1.0u"K") ∘ TVExp(),
            K2 = TVScale(1.0u"K") ∘ TVShift(216.0) ∘ TVScale(45.0) ∘ TVLogistic(),)
# plot(Tran, FulbrightPores.k_f_nt(transform(as(sigma_nt), zeros(4))).(Tran), yscale=:log10, label="MYEGA, def")
trans_myega = as(merge(porous_nt, (;k_f=as(myega_nt))))
trans_sigma = as(merge(porous_nt, (;k_f=as(sigma_nt))))

remap_vars(fitprm) = vcat(repeat(fitprm[1:2], 4), fitprm[3:end])
trans_nm = as((separate=as(Vector, as((N=NTV, S=STV)), 4),
                shared=as((C=CTV, Rs=RsTV, k_f=as(myega_nt))),))
trans_ns_sem_s = as((C=CTV, Rs=RsTV, k_f=as(sigma_nt)))
trans_ns_sem = as((separate=as(Vector, as((C=CTV,)), 4),
                shared=as((Rs=RsTV, k_f=as(sigma_nt))),))
trans_ns = as((separate=as(Vector, as((N=NTV, S=STV)), 4),
                shared=as((C=CTV, Rs=RsTV, k_f=as(sigma_nt))),))


function badNS(nt) 
    if hasproperty(nt, :N) && hasproperty(nt, :S)
        return (nt.N(0) < 1/nt.S(0))
    end
    error("Unhandled nt in badNS")
end
objnf_mc = OptimizationFunction((x,y)->objn_pd(x,y,tweight=1e-1, badprms=badNS), AutoForwardDiff())

@time opt1234myega = solve(OP(objnf_mc, zeros(14), (trans_ns, cases_Rp.pomc, cases_Rp.pdfit)), optalg)

# p0_sem = vcat(repeat([-2.5], 4), zeros(4))
p0_sem = zeros(9)
blankplot_hrC()
modconvtplot!(LyoPronto.gen_nsol_pd(p0_sem, trans_ns_sem, cases_Rp.pomc, cases_Rp.pdfit)..., c=[1 2 3 4])

@time opt1234myegasem = solve(OP(objnf_pd, p0_sem, (trans_ns_sem, cases_Rp.pomc, cases_Rp.pdfit)), optalg, show_trace=true)
@time opt1234myegasems = solve(OP(objnf_pd, zeros(6), (trans_ns_sem_s, cases_Rp.pomc, cases_Rp.pdfit)), optalg, show_trace=true)

fit_mc = transform(trans_ns, opt1234myega.u)
fit_mc_sem = transform(trans_ns_sem, opt1234myegasem.u)
fit_mc_sem_s = transform(trans_ns_sem_s, opt1234myegasems.u)
# all_pomc = [setproperties(pomc, merge(sep, fit_mc.shared)) for (pomc, sep) in zip(cases_Rp.pomc, fit_mc.separate)]
# all_pomc = [setproperties(pomc, merge(sep, fit_mc_sem.shared)) for (pomc, sep) in zip(cases_Rp.pomc, fit_mc_sem.separate)];
all_pomc = [setproperties(pomc, fit_mc_sem_s) for pomc in cases_Rp.pomc];

sols_m = [solve(ODEProblem(pomc), Rodas3()) for pomc in all_pomc];
cases_both = Table(cases_Rp, (pomc=all_pomc,sol_mc=sols_m))

begin
myega_plots = map(cases_both, psols) do case, solp
    @info "check" obj_expT(case.sol_mc, case.pdfit, tweight=0) obj_expT(case.sol_mc, case.pdfit, tweight=1)
    blankplot_hrC(legend=:bottomright, title=case.name)
    plot!(Tsh, label=L"T_{sh}", c=:black)
    modconvtplot!(case.sol_mc)
    modconvtplot!(solp, c=:green, alpha=0.6, label=L"$T_f$, Pikal")
    plot!(case.pdfit)
    plot!(xlim=(-0.5, 21), ylim=(-45, -2))
end
plot(myega_plots..., size=(600, 500))
end
savefig(plotsdir("dextran_T_semC.svg"))
savefig(plotsdir("dextran_T_semC.pdf"))

l = range(0u"cm", hf0, length=100)
begin
blankplot_lRp()
for (i, case) in enumerate(cases_both)
    # plot!(subplot=4, l, case.Rp.(l), ls=:dash, c=parse(Int, case.name), label="3PP $(case.name)")
    mcrpplot!(case.sol_mc, c=i, label="PBM $(case.name)")
    hdT, RpT = case.hRp
    plot!(l, case.po.Rp.(l), ls=:dash, c=i, label="3PP $(case.name)")
    scatter!(hdT, RpT, alpha=0.5, c=i, label="Rp from T, $(case.name)", seriestype=:samplemarkers)
end
plot!(xlim=(0, 0.4), ylim=(0, 15))
plot!()
end
savefig(plotsdir("dextran_Rp.svg"))
savefig(plotsdir("dextran_Rp.pdf"))



ns_SEM = Table(map(cases_Rp, fit_mc.separate, fit_mc_sem.separate) do case, fit, fitsem
    # μ1 = π/6*sum(case.Df.f .* case.Df.D.^3)/sum(case.Df.f)
    # N_sem = 1/μ1 |> u"mL^-1"
    # S_sem = π^2/36 * sum(case.Df.f .* case.Df.D.^6)/sum(case.Df.f) / Vfill |> u"mL"
    N_sem = case.pomc.N(0) |> u"mL^-1"
    S_sem = case.pomc.S(0) |> u"mL"
    C_sem = fitsem.C

    # N_mod = case.sol_mc[3,end] / Vfill
    # S_mod = case.sol_mc[4,end]*u"mm"^6 / Vfill |> u"mL"
    N_fit = fit.N(0) |> u"mL^-1"
    S_fit = fit.S(0) |> u"mL"
    C_fit = fit_mc.shared.C
    C_shared = fit_mc_sem_s.C
    # @info "check" N S sum(case.Df.f)
    # return (;N_sem,N_mod, N_fit, S_sem, S_mod, S_fit, C_sem, C_fit)
    return (;N_sem, N_fit, S_sem, S_fit, C_fit, C_sem, C_shared)
end)

# @df ns_SEM scatter(:N_sem, :S_sem, scale=:log10)
@df ns_SEM scatter((:N_sem)', (:S_sem)', label=permutedims(cases_both.name), scale=:log10)
plot!([10^6.4, 10^7.4], [1e-6, 1e-7], ls=:dash, label="slope -1", c=:gray)
# @df ns_SEM scatter(:N_sem, :C_sem, scale=:log10)
@df ns_SEM scatter((:N_sem)', (:C_sem)', label=permutedims(cases_both.name), scale=:log10)
# plot!([10^6.4, 10^7.4], [10^-0.5, 10^-1.5], label="slope -1")
@df ns_SEM scatter((:S_sem)', (:C_sem)', label=permutedims(cases_both.name), scale=:log10)
plot!(xlabel="S", ylabel="C")
# plot!([1e-15, 1e-13], [10^-1.5, 10^-0.5], label="slope 2")
