
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
    (;name=case, tT, Df = readcaseD.(case))
end)

# Hard-coded for A3, A5, A6, B
t_end = [(18.0u"hr", 24.0u"hr"),
        (19.0u"hr", 26.0u"hr"),
        (20.0u"hr", 26.0u"hr"),
        (17.0u"hr", 26.0u"hr")]
csolid = [0.05u"g/mL", 0.15u"g/mL", 0.2u"g/mL", 0.05u"g/mL"]

cases = Table(data, (;t_end, csolid))

pdends = [13u"hr", 12.5u"hr", 14u"hr", 14u"hr"]
cases_pdfit = Table(map(cases, pdends) do case, pdend
    (; tT, t_end) = case
    pdfit = @df filter(x->x.t<pdend, tT) PrimaryDryFit(:t .- :t[1], :T, t_end .- :t[1])
    merge(case, (;pdfit))
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

Ap, Av = π.*get_vial_radii("10R").^2 .|> u"cm^2"
Vfill = 3u"mL"
hf0 = Vfill/Ap
Tsh = RampedVariable([-40.0, -5]u"°C" .|> u"K", 1u"K/minute")
pch = RampedVariable(10.0u"Pa")
Kshf_g = ConstPhysProp(10.0u"W/m^2/K")
Rp_g = RpFormFit(0.1u"cm^2*Torr*hr/g", 14u"cm*Torr*hr/g", 0.1u"cm^-1")
ρsolution = 1.0u"g/mL"
csolid = 0.05u"g/mL"

po_base = ParamObjPikal(Rp_g, hf0, csolid, ρsolution, Kshf_g, Av, Ap, pch, Tsh) 
sol = solve(ODEProblem(po_base), Rodas3())
blankplot_hrC()
modconvtplot!(sol)
plot!(cases_pdfit.pdfit[1], label="fit", color=:black, linewidth=2)

# Fit to get Kshf
trans_KRp = KRp_transform_basic(Kshf_g.val, Rp_g.R0, Rp_g.A1, Rp_g.A2)

objf_pd = OptimizationFunction((x,y)->obj_pd(x,y,tweight=1e-2), AutoForwardDiff())
@time objf_pd(fill(0.0, 4), (trans_KRp, po_base, cases_pdfit.pdfit[1]))

KRp_opts = map(cases_pdfit) do case
    po = @set po_base.csolid = case.csolid
    opt_KRp = solve(OP(objf_pd, fill(0.0, 4), (trans_KRp, po, case.pdfit)), optalg)
    fitprm = transform(trans_KRp, opt_KRp.u)
    sol = gen_sol_pd(opt_KRp.u, trans_KRp, po)
    # @info "check" fitprm obj_pd(opt_KRp.u, (trans_KRp, po, case.pdfit), verbose=true) sol.t[end]
    blankplot_hrC(legend=:bottomright)
    modconvtplot!(sol)
    plot!(case.pdfit)
end
plot(KRp_opts...)

# Based on fit results... hardcoding this might be dangerous but whatever
@reset po_base.Kshf = ConstPhysProp(21.0u"W/m^2/K")

