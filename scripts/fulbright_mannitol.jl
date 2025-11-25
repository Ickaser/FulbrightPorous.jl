fp_plot_defaults()

begin
Kv = load(datadir("exp_pro", "gravimetric_kv.jld2"))["Kv"]

dat1raw = CSV.read(datadir("exp_raw", "Fulbright2025", "2025-05-21-10_MFD_ISW.csv"), Table; header=8)
dat1pd = filter(x->x.Phase==4, dat1raw)
dat1 = map(x->microfd_fulbright_columnrename(x, dat1pd.CycleTime[1]), dat1pd)
overnightfix(dat1)

# plot(dat1raw.Phase)

# t_end1 = identify_pd_end(dat1.t, dat1.pch_pir, Val(:2der), window_width=91,  tmin=1u"hr",  tmax=50u"hr")
t_end1 = identify_pd_end(dat1.t, dat1.pch_pir, Val(:onoff), window_width=91,  tmin=1u"hr",  tmax=50u"hr")

pch1 = RampedVariable(100.0u"mTorr")
# @df dat1 plot(:t, :pch_cm, ylim=(0, 200))
# @df dat1 plot!(:t, :pch_pir, ylim=(0, 200))
# plot!(pch1)
# tendplot!(t_end1)

Tsh1 = RampedVariable(uconvert.(u"K", [-45.0, -10]u"°C"), 0.5u"K/minute")
# @df dat1 plot(:t, :Tsh_i, xlim=(0, 2))
# plot!(Tsh1)

pl1 = blankplot_hrC()
@df dat1 exptfplot!(:t, :T1, :T2, :T3, :T4)
tendplot!(t_end1)

fitdat1 = @df dat1 PrimaryDryFit(:t[0u"hr".< :t] , 
    # (:T1[0u"hr" .< :t .< 3.7u"hr"],
    #  :T2[0u"hr" .< :t .< 3.7u"hr"],
    #  :T3[0u"hr" .< :t .< 3.7u"hr"],
    #  :T4[0u"hr" .< :t .< 3.7u"hr"],
    # ),
    (:T1[0u"hr" .< :t .< 7u"hr"],
     :T2[0u"hr" .< :t .< 7u"hr"],
     :T3[0u"hr" .< :t .< 7u"hr"],
     :T4[0u"hr" .< :t .< 7u"hr"],
    );
    t_end=t_end1
)
plot(fitdat1)

dat2raw = CSV.read(datadir("exp_raw", "Fulbright2025", "2025-05-23-10_MFD_ISW.csv"), Table; header=8)
dat2pd = filter(x->x.Phase==4, dat2raw)[begin:(end÷5)]
dat2 = map(x->microfd_fulbright_columnrename(x, dat2pd.CycleTime[1]), dat2pd)
overnightfix(dat2)

# plot(dat2raw.Phase)

# t_end2 = identify_pd_end(dat2.t, dat2.pch_pir, Val(:2der), window_width=91,  tmin=1u"hr",  tmax=50u"hr")
t_end2 = identify_pd_end(dat2.t, dat2.pch_pir, Val(:onoff), window_width=91,  tmin=1u"hr",  tmax=50u"hr")

pch2 = RampedVariable(100.0u"mTorr")
# @df dat2 plot(:t, :pch_cm, ylim=(0, 200))
# @df dat2 plot!(:t, :pch_pir, ylim=(0, 200))
# plot!(pch2)
# tendplot!(t_end2)

Tsh2 = RampedVariable(uconvert.(u"K", [-45.0,  10]u"°C"), 0.5*55/35*u"K/minute")
# @df dat2 plot(:t, :Tsh_i, xlim=(0, 2))
# plot!(Tsh2)

pl1 = blankplot_hrC()
# @df dat2 exptfplot!(:t, :T1, :T2, :T3, :T4)
# tendplot!(t_end2)

fitdat2 = @df dat2 PrimaryDryFit(:t[0u"hr".< :t] , 
    # (:T1[0u"hr" .< :t .< 6u"hr"],
    #  :T2[0u"hr" .< :t .< 6.5u"hr"],
    #  :T3[0u"hr" .< :t .< 2u"hr"],
    #  :T4[0u"hr" .< :t .< 6.5u"hr"],
    (:T1[0u"hr" .< :t .< 4u"hr"],
     :T2[0u"hr" .< :t .< 4u"hr"],
     :T3[0u"hr" .< :t .< 4u"hr"],
     :T4[0u"hr" .< :t .< 4u"hr"],
    );
    t_end=t_end2
)
plot(fitdat2)

end

blankplot_hrC()
plot!(fitdat1, c=1)
plot!(fitdat2, c=2)

begin
rad_i, rad_o = get_vial_radii("10R")
A_p = π*rad_i^2
A_v = π*rad_o^2
fillvol = 2.0u"mL"
hf0 = fillvol/A_p |> u"cm"

Rp_g = RpFormFit(1.4u"cm^2*Torr*hr/g", 16.0u"cm*Torr*hr/g", 0.0u"1/cm")

csolid = 0.05u"g/mL"
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
    (Kv, A_v, A_p),
    (pch1, Tsh1),
])
pomc1 = ParamObjMC([
    (hf0, csolid, ρsolution),
    (Kv, A_v, A_p),
    (pch1, Tsh1),
    (ConstPhysProp(k0), C, N, S, Rs)
])

Tsh2 = RampedVariable(uconvert.(u"K", [-45.0, 10]u"°C"), 0.5u"K/minute")
pop2 = @set pop1.Tsh = Tsh2
# @reset pop2.Rp.R0 = 1e-8u"cm^2*Torr*hr/g"
# @reset pop2.Rp.A1 = 32.3u"cm*Torr*hr/g"
pomc2 = @set pomc1.Tsh = Tsh2
end

# blankplot_hrC()
# modconvtplot!(solve(ODEProblem(pop2), LyoPronto.odealg_chunk2))

cases = Table([ (name="1", pop=pop1, pomc=pomc1, fitdat=fitdat1, dat=dat1),
                (name="2", pop=pop2, pomc=pomc2, fitdat=fitdat2, dat=dat2),
])

trans_Rp = Rp_transform_basic(1.4u"cm^2*Torr*hr/g", 16.0u"cm*Torr*hr/g", 1e-2u"1/cm")
obj_Rpf = OptimizationFunction((x,y)->obj_pd(x,y,tweight=1e-1), AutoForwardDiff(chunksize=3))
@time obj_Rpf([0.0, 1.0, 0.0], (trans_Rp, pop1, fitdat1))

single_Rps = map(cases) do case
    (;name, pop, fitdat, dat) = case
    u0 = zeros(3)
    @time opt = solve(OP(obj_Rpf, u0, (trans_Rp, pop, fitdat)), optalg)#, g_tol=1e-6)
    Rp = RpFormFit(transform(trans_Rp, opt.u).Rp...)
    sol = gen_sol_pd(opt.u, trans_Rp, pop, fitdat)
    @info "analyis" name opt Rp
    hRps = map(x->calc_hRp_T(case.pop, case.fitdat, i=x), 1:4)
    return (opt=opt, sol_p=sol, Rp=Rp, hRps=hRps)
end

obj_n = OptimizationFunction((x,y)->LyoPronto.objn_pd(x,y,tweight=1e-1), AutoForwardDiff())
rp12 = @time solve(OP(obj_n, fill(0.0, 3), (trans_Rp, cases.pop[1:2], cases.fitdat[1:2])), optalg;)

cases_Rps = Table(cases, single_Rps)

begin
pls = map(cases_Rps) do case 
    plT = blankplot_hrC()
    # plot!(case.fitdat, title=case.name)
    @df case.dat exptfplot!(:t, :T1, :T2, :T3, :T4, labsuffix=", full",)
    tendplot!(case.fitdat.t[case.fitdat.Tf_iend[1]], label="end useful")
    tendplot!(case.fitdat.t_end)
    modconvtplot!(case.sol_p, label="3PP fit, $(case.name)")
    modconvtplot!(gen_sol_pd(rp12.u, trans_Rp, case.pop), c=:purple, label="3PP fit, 1&2")
    plot!(legend = :bottomright)
    plot!(case.pop.Tsh, tmax=case.fitdat.t_end[end], c=:black, label=L"T_\mathrm{sh}")
    return plT
end
plot(pls..., layout=(1,2), size=(800, 400))
end

l = range(0u"cm", hf0, length=100)
begin
marks = [:o, :x, :square, :diamond]
pls = map(cases_Rps) do case
    plr = blankplot_lRp()
    for i in 1:4
        hdT, RpT = case.hRps[i]
        plot!(hdT, RpT, c=:gray, marker=marks[i], offset=mod1(i, 4)*5, seriestype=:samplemarkers, step=20, linealpha=0.8, label="Rp from T, $i")
    end
    plot!(l, case.Rp.(l), c=1, label="3PP fit directly to T")
    plot!(l, RpFormFit(transform(trans_Rp, rp12.u).Rp...).(l), c=2, label="3PP fit to both T")
    plot!(legend = :bottomright)
    return plr
end
plot(pls..., layout=(1,2), xlim=(0, 0.3), ylim=(0, 20), size=(800, 400))
end
