fp_plot_defaults()
Latexify.set_default(labelformat=:square)

begin
Kv = load(datadir("exp_pro", "gravimetric_kv.jld2"))["Kv"]

MWsuc = 342.3u"g/mol"
MWwat = 18.01u"g/mol"
function freezepoint(csolid)
    xsol = csolid/MWsuc / (csolid/MWsuc + 1u"g/mL"/MWwat)
    ΔHm = 6.01u"kJ/mol"
    Tm_0 = 273.15u"K"
    ΔT = xsol*u"R"*Tm_0^2/ΔHm
    return Tm_0 - ΔT
end

freezepoint(0.05u"g/mL")

dat1raw = CSV.read(datadir("exp_raw", "Fulbright2025", "2025-01-20-10_MFD_ISW.csv"), Table; header=8)
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

Tsh1 = RampedVariable(uconvert.(u"K", [-47.0, -20]u"°C"), 0.5u"K/minute")
# @df dat1 plot(:t, :Tsh_i, xlim=(0, 2))
# plot!(Tsh1)

@df dat1 exptfplot(:t, :T1, :T2, :T3, :T4)

fitdat1 = @df dat1 ExpFitData(:t, 
    TfData(:T1[:t .< 4.8u"hr"]),
    TfData(:T2[:t .< 5u"hr"]),
    TfData(:T3[:t .< 5u"hr"]),
    TfData(:T4[:t .< 5.7u"hr"]),
    EndTimeData(t_end1)
)
plot(fitdat1)

dat2raw = CSV.read(datadir("exp_raw", "Fulbright2025", "2025-01-28-10_MFD_ISW.csv"), Table; header=8)
dat2pd = filter(x->x.Phase==4, dat2raw)
dat2 = map(x->microfd_fulbright_columnrename(x, dat2pd.CycleTime[1]), dat2pd)
overnightfix(dat2)
# begin
# @df filter(x->x.Phase == 1, dat2raw) exptfplot(:CycleTime, :TP1, :TP2, :TP3, :TP4, sampmarks=true)
# plot!(xlim=(Dates.Time("11:30"), Dates.Time("11:55")))
# plot!(ylim=(-10, 0))
# end
# dat2.T1 .+= 1u"K" # Fudge to account for TC calibration
# dat2.T4 .-= 0.5u"K"

# t_end2 = identify_pd_end(dat2.t, dat2.pch_pir, Val(:2der), window_width=91,  tmin=1u"hr",  tmax=50u"hr")
t_end2 = identify_pd_end(dat2.t, dat2.pch_pir, Val(:onoff), window_width=91,  tmin=1u"hr",  tmax=50u"hr")

pch2 = RampedVariable(100.0u"mTorr")
# @df dat2 plot(:t, :pch_cm, ylim=(0, 200))
# @df dat2 plot!(:t, :pch_pir, ylim=(0, 200))
# plot!(pch2)
# tendplot!(t_end2)

Tsh2 = RampedVariable(uconvert.(u"K", [-47.0, -10]u"°C"), 0.5u"K/minute")
# @df dat2 plot(:t, :Tsh_i, xlim=(0, 2))
# plot!(Tsh2)

# @df dat2 exptfplot(:t, :T1, :T2, :T3, :T4)

fitdat2 = @df dat2 ExpFitData(:t,
    TfData(:T1[:t .< 7u"hr"]),
    TfData(:T2[:t .< 7u"hr"]),
    TfData(:T3[:t .< 4u"hr"]),
    TfData(:T4[:t .< 3u"hr"]),
    EndTimeData(t_end2)
)
# plot(fitdat2)

dat3raw = CSV.read(datadir("exp_raw", "Fulbright2025", "2025-02-03-10_MFD_ISW.csv"), Table; header=8)
dat3pd = filter(x->x.Phase==4, dat3raw)
dat3 = map(x->microfd_fulbright_columnrename(x, dat3pd.CycleTime[1]), dat3pd)
overnightfix(dat3)

# t_end3 = identify_pd_end(dat3.t, dat3.pch_pir, Val(:2der), window_width=91,  tmin=1u"hr",  tmax=50u"hr")
t_end3 = identify_pd_end(dat3.t, dat3.pch_pir, Val(:onoff), window_width=91,  tmin=1u"hr",  tmax=50u"hr")

pch3 = RampedVariable(100.0u"mTorr")
# @df dat3 plot(:t, :pch_cm, ylim=(0, 200))
# @df dat3 plot!(:t, :pch_pir, ylim=(0, 200))
# plot!(pch3)
# tendplot!(t_end3)

Tsh3 = RampedVariable(uconvert.(u"K", [-47.0, 0.0]u"°C"), 0.5u"K/minute")
# @df dat3 plot(:t, :Tsh_i, xlim=(0, 2))
# plot!(Tsh3)

# @df dat3 exptfplot(:t, :T1, :T2, :T3, :T4)

fitdat3 = @df dat3 ExpFitData(:t,
    TfData(:T1[:t .< 3u"hr"]),
    TfData(:T2[:t .< 4u"hr"]),
    TfData(:T3[:t .< 3.5u"hr"]),
    TfData(:T4[:t .< 3.5u"hr"]),
    EndTimeData(t_end3)
)
# plot(fitdat3)


# ---------------
# 10% Sucrose cycles 
dat4raw = CSV.read(datadir("exp_raw", "Fulbright2025", "2025-02-25-10_MFD_IW.csv"), Table; header=8)
dat4pd = filter(x->x.Phase==4, dat4raw)
dat4 = map(x->microfd_fulbright_columnrename(x, dat4pd.CycleTime[1]), dat4pd)
overnightfix(dat4)

# t_end4 = identify_pd_end(dat4.t, dat4.pch_pir, Val(:2der), window_width=91,  tmin=1u"hr",  tmax=50u"hr")
t_end4 = identify_pd_end(dat4.t, dat4.pch_pir, Val(:onoff), window_width=91,  tmin=1u"hr",  tmax=50u"hr")

pch4 = RampedVariable(100.0u"mTorr")
# @df dat4 plot(:t, :pch_cm, ylim=(0, 200))
# @df dat4 plot!(:t, :pch_pir, ylim=(0, 200))
# plot!(pch4)
# tendplot!(t_end4)

Tsh4 = RampedVariable(uconvert.(u"K", [-47.0, -20.0]u"°C"), 0.5u"K/minute")
# @df dat4 plot(:t, :Tsh_i, xlim=(0, 2))
# plot!(Tsh4)

# @df dat4 exptfplot(:t, :T1, :T2, :T3, :T4)

fitdat4 = @df dat4 ExpFitData(:t,
    TfData(:T1[:t .< 6.4u"hr"]),
    TfData(:T2[:t .< 6.5u"hr"]),
    TfData(:T3[:t .< 7u"hr"]),
    TfData(:T4[:t .< 8.5u"hr"]),
    EndTimeData(t_end4)
)
plot(fitdat4)

dat5raw = CSV.read(datadir("exp_raw", "Fulbright2025", "2025-02-27-10_MFD_IW.csv"), Table; header=8)
dat5pd = filter(x->x.Phase==4, dat5raw)
dat5 = map(x->microfd_fulbright_columnrename(x, dat5pd.CycleTime[1]), dat5pd)
overnightfix(dat5)

# t_end5 = identify_pd_end(dat5.t, dat5.pch_pir, Val(:2der), window_width=91,  tmin=1u"hr",  tmax=50u"hr")
t_end5 = identify_pd_end(dat5.t, dat5.pch_pir, Val(:onoff), window_width=91,  tmin=1u"hr",  tmax=50u"hr")

pch5 = RampedVariable(100.0u"mTorr")
# @df dat5 plot(:t, :pch_cm, ylim=(0, 200))
# @df dat5 plot!(:t, :pch_pir, ylim=(0, 200))
# plot!(pch5)
# tendplot!(t_end5)

Tsh5 = RampedVariable(uconvert.(u"K", [-47.0, -30.0]u"°C"), 0.5u"K/minute")
# @df dat5 plot(:t, :Tsh_i, xlim=(0, 2))
# plot!(Tsh5)

# @df dat5 exptfplot(:t, :T1, :T2, :T3, :T4)

fitdat5 = @df dat5 ExpFitData(:t,
    TfData(:T1[:t .< 11u"hr"] .-2u"K"),
    TfData(:T2[:t .< 14u"hr"] .-2u"K"),
    TfData(:T3[:t .< 13u"hr"] .-2u"K"),
    TfData(:T4[:t .< 10u"hr"] .-2u"K"),
    EndTimeData(t_end5)
)
# plot(fitdat5)

dat6raw = CSV.read(datadir("exp_raw", "Fulbright2025", "2025-05-14-10_MFD_ISW.csv"), Table; header=8)
dat6pd = filter(x->x.Phase==4, dat6raw)
dat6 = map(x->microfd_fulbright_columnrename(x, dat6pd.CycleTime[1]), dat6pd)
overnightfix(dat6)

# t_end5 = identify_pd_end(dat5.t, dat5.pch_pir, Val(:2der), window_width=91,  tmin=1u"hr",  tmax=50u"hr")
t_end6 = identify_pd_end(dat6.t, dat6.pch_pir, Val(:onoff), window_width=91,  tmin=1u"hr",  tmax=50u"hr")

pch6 = RampedVariable(100.0u"mTorr")
# @df dat6 plot(:t, :pch_cm, ylim=(0, 200))
# @df dat6 plot!(:t, :pch_pir, ylim=(0, 200))
# plot!(pch6)
# tendplot!(t_end6)

Tsh6 = RampedVariable(uconvert.(u"K", [-45.0, -0]u"°C"), 0.5u"K/minute")
# @df dat6 plot(:t, :Tsh_i, xlim=(0, 2))
# plot!(Tsh6)

@df dat6 exptfplot(:t, :T1, :T2, :T3, :T4)

fitdat6 = @df dat6 ExpFitData(:t,
    TfData(:T1[:t .< 3u"hr"]),
    TfData(:T2[:t .< 3u"hr"]),
    TfData(:T3[:t .< 2.5u"hr"]),
    TfData(:T4[:t .< 3u"hr"]),
    EndTimeData(t_end6)
)
# plot(fitdat6)

dat7raw = CSV.read(datadir("exp_raw", "Fulbright2025", "2025-05-20-10_MFD_ISW.csv"), Table; header=8)
dat7pd = filter(x->x.Phase==4, dat7raw)
dat7 = map(x->microfd_fulbright_columnrename(x, dat7pd.CycleTime[1]), dat7pd)
overnightfix(dat7)

# t_end5 = identify_pd_end(dat5.t, dat5.pch_pir, Val(:2der), window_width=91,  tmin=1u"hr",  tmax=50u"hr")
t_end7 = identify_pd_end(dat7.t, dat7.pch_pir, Val(:onoff), window_width=91,  tmin=1u"hr",  tmax=50u"hr")

pch7 = RampedVariable(100.0u"mTorr")
# @df dat7 plot(:t, :pch_cm, ylim=(0, 200))
# @df dat7 plot!(:t, :pch_pir, ylim=(0, 200))
# plot!(pch7)
# tendplot!(t_end7)

Tsh7 = RampedVariable(uconvert.(u"K", [-46.0, -20]u"°C"), 0.5u"K/minute")
# @df dat7 plot(:t, :Tsh_i, xlim=(0, 2))
# plot!(Tsh7)

# @df dat7 exptfplot(:t, :T1, :T2, :T3, :T4)
@df dat7raw plot([:TP1, :TP2, :TP3, :TP4], xlim=(0, 120), ylim=(-5, 0))
# # @df dat7raw plot(:TP1 .- :TP2)

fitdat7 = @df dat7 ExpFitData(:t,
    TfData(:T1[:t .< 6.5u"hr"]),
    TfData(:T2[:t .< 7.5u"hr"]),
    TfData(:T3[:t .< 9u"hr"]),
    TfData(:T4[:t .< 8u"hr"]),
    EndTimeData(t_end7)
)
# plot(fitdat7)

end;




blankplot_hrC()
plot!(fitdat1, c=1, nmarks=15)
plot!(fitdat2, c=2, nmarks=15)
plot!(fitdat3, c=3, nmarks=15)
# plot!(fitdat6, c=6)
plot!(fitdat7, c=7, nmarks=15)
blankplot_hrC()
plot!(fitdat4, c=4, nmarks=15)
plot!(fitdat5, c=5, nmarks=15)
# plot!(xlim=(0, 10), ylim=(-36,-30))

# ------------- Check Pirani

# @df dat1raw plot(:t, :VacPirani, ylim=(90, 200))
# @df dat2raw plot!(:VacPirani, )
# @df dat3raw plot!(:VacPirani, )
# @df dat4raw plot!(:VacPirani, )
# @df dat5raw plot!(:VacPirani, ylim=(90, 180))



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

Rp_g = RpFormFit(0.8u"cm^2*Torr*hr/g", 14.0u"cm*Torr*hr/g", 1.0u"1/cm")

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

pop2 = @set pop1.Tsh = Tsh2
pop3 = @set pop1.Tsh = Tsh3
pop7 = @set pop1.Tsh = Tsh7
pomc2 = @set pomc1.Tsh = Tsh2
pomc3 = @set pomc1.Tsh = Tsh3
pomc7 = @set pomc1.Tsh = Tsh7

pop4 = @set pop1.csolid = 0.1u"g/mL"
pop5 = @set pop4.Tsh = Tsh5
pomc4 = @set pomc1.csolid = 0.1u"g/mL"
pomc5 = @set pomc4.Tsh = Tsh5
end


cases = Table([ (name="1", pop=pop1, pomc=pomc1, fitdat=fitdat1, dat=dat1),
                (name="2", pop=pop2, pomc=pomc2, fitdat=fitdat2, dat=dat2),
                (name="3", pop=pop3, pomc=pomc3, fitdat=fitdat3, dat=dat3),
                (name="4", pop=pop4, pomc=pomc4, fitdat=fitdat4, dat=dat4),
                (name="5", pop=pop5, pomc=pomc5, fitdat=fitdat5, dat=dat5),
                (name="7", pop=pop7, pomc=pomc7, fitdat=fitdat7, dat=dat7),
                ])

trans_Rp = Rp_transform_basic(0.01u"cm^2*Torr*hr/g", 14u"cm*Torr*hr/g", 1e0u"1/cm")
obj_Rpf = OptimizationFunction((x,y)->obj_pd(x,y,tweight=1e-1, badprms=x->x.Rp.R0<1e-4u"cm^2*Torr*hr/g"), AutoForwardDiff(chunksize=3))
# obj_Rpf = OptimizationFunction(obj_pd, AutoForwardDiff(chunksize=3))

obj_npf = OptimizationFunction((x,y)->LyoPronto.objn_pd(x,y,tweight=1e-1, badprms=x->x.Rp.R0<1e-4u"cm^2*Torr*hr/g"), AutoForwardDiff())

rp123 = @time solve(OP(obj_npf, fill(0.0, 3), (trans_Rp, cases.pop[1:3], cases.fitdat[1:3])), optalg;)


single_Rps = map(cases) do case
    (;name, pop, fitdat, dat) = case
    u0 = if name ∈ ("5","7") 
        [1.0, 0, 1.5] 
    elseif name == "2"
        [1.0, 1.0, 2.0]
    else 
        [1.0, 1.0, 1.0]
    end
    @time opt = solve(OP(obj_Rpf, u0, (trans_Rp, pop, fitdat)), optalg)#, g_tol=1e-6)
    # nls_f = NonlinearFunction(fitdat; tweight=1e-2)
    # @time opt = solve(NonlinearLeastSquaresProblem(nls_f, u0, p=(trans_Rp, pop, fitdat)), LevenbergMarquardt(), reltol=1e-6)
    !SciMLBase.successful_retcode(opt.retcode) && @warn "Unsuccessful fit!!" name 
    Rp = transform(trans_Rp, opt.u).Rp
    sol = gen_sol_pd(opt.u, trans_Rp, pop, fitdat)
    @info "analyis" name opt Rp sol.retcode
    hRps = map(x->calc_hRp_T(case.pop, case.fitdat, i=x), 1:4)
    return (opt=opt.u, sol_p=sol, Rp=Rp, hRps=hRps)
end


cases_Rps = Table(cases, single_Rps)

begin
pls = map(cases_Rps) do case 
    plT = blankplot_hrC()
    plot!(case.fitdat, title=case.name, nmarks=15)
    modconvtplot!(case.sol_p)
    plot!(legend = :bottomright)
    # plot!(case.pop.Tsh, tmax=case.fitdat.t_end, c=:black, label=L"T_\mathrm{sh}")
    return plT
end
plot(pls..., layout=(3,2), size=(800, 1000))
end

begin
l = range(0u"cm", hf0, length=100)
pls = map(cases_Rps) do case
    plr = blankplot_lRp()
    for i in 1:4
        hdT, RpT = case.hRps[i]
        plot!(hdT, RpT, c=i, label="Rp from T, $i")
    end
    plot!(l, case.Rp.(l), label="3PP fit directly to T")
    plot!(legend = :bottomright)
    return plr
end
plot(pls..., layout=(3,2), size=(800, 1000))
end

function make_mcpass(casenums)
    subdat = filter(x->parse(Int, x.name)∈casenums, cases)
    (subdat.pomc, subdat.fitdat)
end

pass123 = make_mcpass([1, 2, 3])
pass45 = make_mcpass([4, 5])
pass7 = make_mcpass([7])

porous_nt = FulbrightPorous.porous_transform_bounded(fillvol)
# myega_nt = FulbrightPorous.myega_transform_basic(1e42u"s^-1", 2.0e4u"K", 1e2u"K")
sigma_nt = FulbrightPorous.sigmoid_transform_bounded(1e-10u"s^-1", 1.0u"K")
# trans_myega = merge(porous_nt, myega_nt)
trans_sigma = merge(porous_nt, sigma_nt)

# objn_myega(fitprm, pf) = objn_mc(fitprm, (trans_myega, pf...), tweight=1e-2)
obj_mc(fitprm, tpf) = objn_pd(fitprm, tpf; badprms=FulbrightPorous.badNS, tweight=1e-2)
objf_mc = OptimizationFunction(obj_mc, AutoForwardDiff(chunksize=4))
objf_mc(rand(7), (trans_myega, pass123...))

p0_sigma = [-1, 0, -1, 0.1, 0, 3, -1.5, 1.0]

## Alternate analyses
@time opt123sigm = solve(OP(objf_mc, p0_sigma, (trans_sigma, pass123...)), optalg)
@time opt45sigm = solve(OP(objf_mc, p0_sigma, (trans_sigma, pass45...)), optalg)
@time opt7sigm = solve(OP(objf_mc, p0_sigma, (trans_sigma, pass7...)), optalg, time_limit=15)

# p0_myega = [0, 0, -1, 0.1, 1, 0, 0]
# @time opt123myega = solve(OP(objf_mc, p0_myega, (trans_myega, pass123...)), optalg, time_limit=15, g_tol=1e-6)
# @time opt45myega = solve(OP(objf_mc, p0_myega, (trans_myega, pass45...)), optalg, time_limit=15, g_tol=1e-6)

## Desired analysis

remap_vars(fitprm) = vcat(fitprm[1:2], fitprm[1:2], fitprm[1:2], fitprm[3:end])
trans_2s = as((separate=as(Vector, as((N=porous_nt.N, S=porous_nt.S)), 3),
                sep_inds=Constant([1, 1, 1, 2, 2, 3]),
                shared=merge(as((C=porous_nt.C, Rs=porous_nt.Rs)), sigma_nt),))
# trans_2m = @set trans_2s.shared.k_f = myega_nt.k_f
# trans_2s_k = @set trans_2s.shared = merge(trans_2s.shared, K_transform_basic(Kv.val))
function gen_pbmsols(opt, tr)
    sols = if hasproperty(tr, :shared)
        gen_nsol_pd(opt, tr, cases.pomc)
    else
        map(cases_Rps) do case
            sol = gen_sol_pd(opt, tr, case.pomc)
        end
    end
    return Table((sol_mc=sols,))
end


# p0_2m = remap_vars(p0_myega)
# p0_2m = remap_vars(opt123myega)
# @time opt_2m = solve(OP(objf_mc_c5, p0_2m, (trans_2m, cases_Rps.pomc, cases_Rps.fitdat)), optalg, time_limit=15)

p0_2s = remap_vars(p0_sigma .+ [-2.0, 1.0, 0, 0, 0, 0, -1, -1])
# p0_2s = remap_vars(opt123sigm.u)
@time opt_2s = solve(OP(objf_mc, p0_2s, (trans_2s, cases_Rps.pomc, cases_Rps.fitdat)), optalg)
# @time opt_2s_k = solve(OP(objf_mc, vcat(opt_2s.u, 0.1), (trans_2s_k, cases_Rps.pomc, cases_Rps.fitdat)), optalg, time_limit=30)


# allopts = [opt123sigm, opt123myega, opt45sigm, opt45myega, opt_2m, opt_2s, opt7sigm, opt_2s_k]
# # check for successes
# [(s = SciMLBase.successful_retcode(opt); !s && @warn "Unsuccessful opt! Check details"; s) for opt in allopts]

# fits = Table((opt=[opt123sigm.u, opt123myega.u, opt45sigm.u, opt45myega.u, opt_2m.u, opt_2s.u, opt7sigm.u, opt_2s_k.u], 
#               tr=[trans_sigma, trans_myega, trans_sigma, trans_myega, trans_2m, trans_2s, trans_sigma, trans_2s_k],
#               cases=[(1,2,3), (1,2,3), (4,5), (4,5), (1,2,3,4,5, 6), (1,2,3,4,5,6), (6,), (1,2,3,4,5,6)],))
fits = Table((opt=[opt123sigm.u, opt45sigm.u, opt7sigm.u, opt_2s.u], 
              tr=[trans_sigma, trans_sigma, trans_sigma, trans_2s],
              cases=[(1,2,3), (4,5), (6,), (1,2,3,4,5, 6)],))
fits = Table(fits, (kf= map(fits) do x
    fitprm = transform(x.tr, x.opt)
    if hasproperty(fitprm, :shared)
        fitprm.shared.k_f
    else
        fitprm.k_f
    end
end ,))
fits = Table(fits, (k_type = map(x->(x.kf isa Sigmoid_k ? Sigmoid_k : MYEGA_k), fits)))

safesave(datadir("exp_pro", "fulbright_sucrose_fitresults.jld2"), Dict("fits"=>fits))

# new_sols = vcat(gen_pbmsols(opt123myega, trans_myega)[1:3], gen_pbmsols(opt45myega, trans_myega)[4:5])
# new_sols = vcat(gen_pbmsols(opt123sigm, trans_sigma)[1:3], gen_pbmsols(opt45sigm, trans_sigma)[4:5], gen_pbmsols(opt7sigm, trans_sigma)[6:6])
# new_sols = gen_pbmsols(opt_2m, trans_2m)
new_sols = gen_pbmsols(opt_2s, trans_2s)
# new_sols = gen_pbmsols(opt_2s_k, trans_2s_k)

cases_both = Table(cases, single_Rps, new_sols)

errs = Table(map(cases_both) do case
    PBM = case.sol_mc
    PP = case.sol_p
    errf = (model, tweight)->obj_expT(model, case.fitdat; tweight)
    PBM_T = sqrt(errf(PBM, 0))*u"K"
    PP_T = sqrt(errf(PP, 0))*u"K"
    t_end_est = ustrip(u"hr", sum(case.fitdat.t_end)/2)
    PBM_t = sqrt(errf(PBM, 1) - errf(PBM, 0))/t_end_est
    PP_t = sqrt(errf(PP, 1) - errf(PP, 0))/t_end_est
    @info "check" case.name PBM_t PBM.t[end] case.fitdat.t_end
    return (;name=case.name, PBM_T, PP_T, PBM_t, PP_t)
end)

# Kv fit: report
# Kv_fit = transform(trans_2s_k, opt_2s_k).shared.Kshf
# @show Kv_fit Kv

# # Plot N and S for different fits
# begin
# i = 1
# # myega_or_sig = [1, 2, 1, 2, 2, 2, 2, 1, 1, 1]
# cs = cgrad(:viridis)
# markers = [:circle, :square, :diamond, :utriangle, :dtriangle, :ltriangle, :rtriangle, :hexagon, :pentagon, :heptagon, :octagon]
# styles = Table((marker=markers, markersize=fill(10, 11)))
# plot(u"μm^-3", u"μm^3", xlabel="N", ylabel="S", xscale=:log10, yscale=:log10, unitformat=latexify)
# pl_fitNS = map(fits[1:3]) do fit # just sigmoid cases
#     fitprm = transform(fit.tr, fit.opt)
#     if hasproperty(fitprm, :shared)
#         N_vals = fitprm.separate .|> x->x.N.val .|> u"μm^-3"
#         S_vals = fitprm.separate .|> x->x.S.val .|> u"μm^3"
#         for (N, S, nums) in zip(N_vals, S_vals, [(1,2,3), (4,5), (7)])
#             scatter!([N], [S], label="$nums, fit to all"; styles[i]...)
#             @info "fits" fitprm.shared.C fit.cases
#             global i += 1
#         end
#     else
#         scatter!([fitprm.N.val], [fitprm.S.val], label="$(fit.cases)"; styles[i]...)
#         @info "fits" fitprm.C fit.cases
#         global i += 1
#     end
# end
# plot!(legend=:right)
# end

begin
pls_T1 = map(cases_both[1:3]) do case
    pl = blankplot_hrC(legend=:bottomright, title="$(case.name), Tsh=$(case.pop.Tsh(Inf*u"s")|>u"°C")")
    plot!(case.fitdat, nmarks=15)
    modconvtplot!(case.sol_mc, c=2, label="PBM: 1&2&3")
    modconvtplot!(case.sol_p, c=:green, linealpha=0.7, label="3PP $(case.name)")
    modconvtplot!(gen_sol_pd(rp123.u, trans_Rp, case.pop), c=:purple, label="3PP 1&2&3")
    plot!(xlim=(0, 16), ylim=(-40, -30), left_margin=20Plots.px)
    plot!(yticks=(-40:2:-30))
    case.name != "1" && plot!(ylabel="", yformatter=:none, left_margin=0Plots.px)
    return pl
end
# plot!(pls_T1[1], ylabel="Temperature [°C]")
plot(pls_T1..., layout=(1,3), size=(700, 400), bottom_margin=20Plots.px)
end
# savefig(plotsdir("fulbright_5sucrose_T.svg"))
# savefig(plotsdir("fulbright_5sucrose_T.pdf"))

begin
errnames = ["5%\n-20", "5%\n-10", "5%\n0", "10%\n-20", "10%\n-30", "5% C\n-20"]
err_permute = [6,1,2,3,4,5]
e1 = @df errs[err_permute] groupedbar(errnames[err_permute], [:PBM_T :PP_T], label=["PBM" "Pikal fit"], ylabel=L"$T_\mathrm{f}$ error", c=[2 3], lw=1, ylim=(0,3), unitformat=:square)
e2 = @df errs[err_permute] groupedbar(errnames[err_permute], [:PBM_t :PP_t], label=["PBM" "Pikal fit"], ylabel=L"$t$ rel. error", c=[2 3], lw=1, ylim=(0,0.15))
# Mark that 10%@-30 is way above plot limit
big_err = string(round(errs[5].PBM_t, digits=2))
annotate!(e2, 5.25, 0.15, "↑ $big_err", font(14, :left, "Computer Modern"))
pls_T1 = map(cases_both[1:3]) do case
    pl = blankplot_hrC(legend=:outerright, ylabel="T", xlabel="Time")
    plot!(case.fitdat, nmarks=10, markeralpha=0.8)
    Tsh = round(typeof(1u"degC"), case.pop.Tsh(Inf*u"s"))
    plot!(title="5%, \n " * L"T_\mathrm{sh} = " * replace(string(Tsh), " "=>""))
    modconvtplot!(case.sol_mc, case.sol_p, c=[2 3], linealpha=0.7, labels=["PBM" "Pikal fit"])
    plot!(xlim=(0, 15.5), ylim=(-40, -30))
    return pl
end
pls_T2 = map(cases_both[4:5]) do case
    pl = blankplot_hrC(legend=:outerright, ylabel="T", xlabel="Time")
    plot!(case.fitdat, nmarks=15, markeralpha=0.8)
    Tsh = round(typeof(1u"degC"), case.pop.Tsh(Inf*u"s"))
    plot!(title="10%, \n " * L"T_\mathrm{sh} = " * replace(string(Tsh), " "=>""))
    modconvtplot!(case.sol_mc, case.sol_p, c=[2 3], linealpha=0.7,labels=["PBM" "Pikal fit"])
    # plot!(xlim=(0, 65), ylim=(-40, -30))
    plot!(ylim=(-40,-30))
    return pl
end
pl_T3 = map(cases_both[6:6]) do case
    pl = blankplot_hrC(legend=:outerright, ylabel="T", xlabel="Time")
    Tsh = round(typeof(1u"degC"), case.pop.Tsh(Inf*u"s"))
    plot!(title="5%, CIN, \n " * L"T_\mathrm{sh} = " * replace(string(Tsh), " "=>""))
    plot!(case.fitdat, nmarks=15)
    modconvtplot!(case.sol_mc, case.sol_p, c=[2 3], linealpha=0.7, labels=["PBM" "Pikal fit"])
    # plot!(xlim=(0, 65), ylim=(-40, -30))
    plot!(ylim=(-40,-30))
    plot!(xlim=(0, 15.5))
    return pl
end
plot(pl_T3..., pls_T1..., pls_T2..., e1, e2, layout=@layout([grid(1,6){0.55h}; [_ grid(1,2){0.8w} _]]), size=(1200, 500))
plot!(legend=:none, top_margin=20Plots.px, bottom_margin=20Plots.px)
for i in [3,4,6]; plot!(subplot=i, ylabel="", yformatter=:none, left_margin=-15Plots.px); end
for i in [2,5]; plot!(subplot=i, ylabel="",left_margin=-5Plots.px); end
plot!(subplot=1, legend=(-0.1,-0.6), legendfontsize=12, left_margin=30Plots.px)
plot!(subplot=8, legend=(1.1,0.5), legendfontsize=12)
end
savefig(plotsdir("fulbright_allsucrose_T.svg"))
# Manually adjust the legend before exporting to PDF

begin
plot(u"°C", u"s^-1", xlabel="Temperature", ylabel="Rate", yscale=:log10, unitformat=:square)
Tran = (range(-40u"°C", -30u"°C", step=0.1u"K")) .|> u"K"
for fit in fits
    # @info "check" fit.kf.(Tran) fit.cases fit.kf
    fit.kf isa Sigmoid_k && plot!(Tran, fit.kf.(Tran), label="$(fit.k_type), $(fit.cases)")
end
# plot!(Tran, transform(trans_myega, p0_myega).k_f.(Tran), label="MYEGA, guess")
# plot!(Tran, transform(trans_sigma, p0_sigma).k_f.(Tran), label="Sigmoid, guess")
# plot!(Tran, transform(myega_nt, zeros(3)).k_f.(Tran), label="MYEGA, def")
# plot!(Tran, transform(sigma_nt, zeros(4)).k_f.(Tran), label="Sigmoid, def")
# plot!(ylim=(1e-16, 1e-5), legend=:bottomright)
plot!(legend=:bottomright)
end

l = range(0u"cm", hf0, length=100)
begin
listed_cases = [1,2,3]
kind = :stats
cs = palette(:Oranges_7)[[3,5,7]]
markercs = palette(:Oranges_7)[[1,4,7]]
markers = [:circle :square :utriangle]
plot(layout=@layout([grid(3,1) a{0.65w}]), size=(800,400))
for (i, case) in enumerate(cases_both[listed_cases])
    Tsh = round(typeof(1u"degC"), case.pop.Tsh(Inf*u"s"))
    for j in 1:4
        hdT, RpT = case.hRps[j]
        nevery = max(1, length(hdT)÷10)
        offset = (j-1)*nevery÷4 + 1
        plot!(subplot = 4, hdT[offset:nevery:end], RpT[offset:nevery:end], alpha=0.5, 
            c=markercs[i], 
            seriestype=:scatter, marker=markers[i],
            label=(j==1 ? L"T_\mathrm{f}(t), T_\mathrm{sh}="*" $(Tsh)" : ""),
            )
    end
end
for (i, ci) in enumerate(listed_cases)
    Tsh = round(typeof(1u"degC"), cases_both[ci].pop.Tsh(Inf*u"s"))
    poreevolplot!(cases_both[ci].sol_mc, label="PBM, $Tsh"; 
        c=cs[i], kind,
        xlims=(0, Inf),
        # seriestype=[:samplemarkers :samplemarkers :samplemarkers :line]
        # marker=[markers[i] markers[i] markers[i] :none], 
        seriestype=:line
    )
end
for (i, ci) in enumerate(listed_cases)
    Tsh = round(typeof(1u"degC"), cases_both[ci].pop.Tsh(Inf*u"s"))
    plot!(subplot=4, l, cases_both[ci].Rp.(l), ls=:dash, c=palette(:Greens_7)[[3 5 7]][i], label="Pikal, $(Tsh)")
end
# poreevolplot!(cases_both[2].sol_mc, label="PBM 2, 1&2&3"; marker=:square, kind)
# poreevolplot!(cases_both[3].sol_mc, label="PBM 3, 1&2&3"; marker=:d, kind)
# poreevolplot!(cases_both[6].sol_mc, label="PBM 7, 7"; marker=:utriangle, kind)
# plot!(subplot=1, xlabel="", yticks=(0:1e9:3e9, string.(0:3, " " ) .* L"\times10^9"), ylim=(0,3e9), left_margin=20Plots.px, bottom_margin=-10Plots.px)
plot!(subplot=1, xlabel="", yscale=:log10, ylim=(5e6, 3e9))
plot!(subplot=2, xlabel="", ylims=(0, Inf), left_margin=35Plots.px, ywiden=1.1, bottom_margin=-10Plots.px, )
plot!(subplot=3, ylims=(0, Inf), ywiden=1.1, bottom_margin=20Plots.px)
plot!(subplot=4, yunitformat=latexify, ylabel="R_p", legend=:outerright, legendfontsize=10, bottom_margin=20Plots.px)
plot!(subplot=4, legend=:outerright)
end
savefig(plotsdir("fulbright_suc_moments.svg"))
savefig(plotsdir("fulbright_suc_moments.pdf"))

# begin
# blankplot_lRp()
# for case in cases_both[[1, 2, 3]]
#     # plot!(subplot=4, l, case.Rp.(l), ls=:dash, c=parse(Int, case.name), label="3PP $(case.name)")
#     mcrpplot!(case.sol_mc, c=parse(Int, case.name), label="PBM mixed-all")
#     for i in 1:4
#         hdT, RpT = case.hRps[i]
#         plot!(hdT, RpT, alpha=0.4, c=parse(Int, case.name), label=(i==1 ? "Rp from T, $(case.name)" : ""))
#     end
#     plot!(l, case.Rp.(l), ls=:dash, c=parse(Int, case.name), label="3PP $(case.name)")
# end
# # plot!(xlim=(0, 0.3), ylim=(0.5, 2.5))
# plot!(size=(400, 400), legend=:outerright)
# end
# savefig(plotsdir("fulbright_sucrose_Rp.svg"))
# savefig(plotsdir("fulbright_sucrose_Rp.pdf"))

# timedistplot(cases_both[5].sol_mc; ywiden=false, scale=:log10, disttype=:beta, legend=:outerright)
begin
plot(u"μm^3", NoUnits, ywiden=false, ylabel="Approx. Prob. Dens.", xlabel="Pore Size", unitformat=latexify)
for case in cases_both[[1, 2, 3, 6, 4, 5]]
    distplot!(case.sol_mc;  disttype=:lognormal, label=case.name)
    plot!(yscale=:log10)
    plot!(xscale=:log10)
end
# plot!(scale=:log10)
plot!()
end
savefig(plotsdir("fulbright_sucrose_poredist.svg"))
savefig(plotsdir("fulbright_sucrose_poredist.pdf"))

res_fit = cases_both[1].sol_mc.prob.p
@info "Best fit params, 1&2&3:" res_fit.N res_fit.S res_fit.C res_fit.Rs res_fit.k_f
res_fit = cases_both[4].sol_mc.prob.p
@info "Best fit params, 4&5:" res_fit.N res_fit.S res_fit.C res_fit.Rs res_fit.k_f
res_fit = cases_both[6].sol_mc.prob.p
@info "Best fit params, 7:" res_fit.N res_fit.S res_fit.C res_fit.Rs res_fit.k_f


## MCMC



# function eval_about_myega(new, opt, pass)
#     objn_myega(new .* opt.u, pass)
# end
# function eval_many(ranges...)
#     ndims = length(ranges)
#     eval_points = [[u..., fill(1.0, 7-ndims)...] for u in Iterators.product(ranges...)]
#     eval_about_myega.(eval_points, [opt123myega], [pass123])
# end

# # spread = [[exp(i), exp(j), fill(1.0,5)...] for i in -5:5, j in -5:5]

# integ_dims(dims...) = dropdims(sum(probs, dims=dims), dims=dims)
# # @time sqerrs = eval_about_myega.(spread, [opt123myega], [pass123])
# begin
# Nr = range(-30, 40, length=11)
# Sr = range(-16, 16, length=12)
# Cr = range(-10, 20, length=5)
# Rsr = range(-2, 10, length=4)
# @time sqerrs = eval_many(Nr, Sr, Cr, Rsr)
# minsqerr = minimum(filter(!isnan, sqerrs)) 
# opt123myega.minimum
# probs = exp.(-map(x->(isnan(x) ? Inf : x-minsqerr), sqerrs))
# probs ./= sum(probs)
# end;
# # contourf(sqerrs[:,:,1,1], cmap=:thermal, lw=0)

# ## Bayesian, Turing
# using Turing

# @model function myega_posterior(opt, pass)

#     Nfac ~ LogUniform(1e-15, 1e20)
#     Sfac ~ LogUniform(1e-6, 1e6)
#     Cfac ~ LogUniform(1e-5, 1e10)
#     Rfac ~ LogUniform(1e-1, 1e5)
#     kfac ~ LogUniform(1e-3, 1e3)
#     Kfac ~ LogUniform(1e-3, 1e3)
#     θfac ~ LogUniform(1e-3, 1e3)

#     prmfacs = [Nfac, Sfac, Cfac, Rfac, kfac, Kfac, θfac]
#     sqerr = eval_about_myega(prmfacs, opt, pass)

#     if !isnan(sqerr)
#         Turing.@addlogprob! -sqerr
#     else
#     #     # @info "check" newpo.N.val newpo.S.val newpo.C.val newpo.Rs.val sol.t[end] sol[1,end]
#         Turing.@addlogprob! -Inf
#     end
# end

# case123_model = myega_posterior(opt123myega, pass123);

# nthreads = 1
# init_set = [exp.(rand(7).-0.5) for i in 1:nthreads]
# @time smallsamp = sample(case123_model, NUTS(), MCMCThreads(), 50, nthreads, initial_params=init_set)
# nthreads = 4
# init_set = [exp.(rand(7).-0.5) for i in 1:nthreads]
# samp_alg = NUTS(100, 0.65; init_ϵ=0.01)
# @time samp = sample(case123_model, samp_alg, MCMCThreads(), 500, nthreads, initial_params=init_set)

# safesave(datadir("exp_pro", "fulbright123_myega_sampling.jld2"), @strdict smallsamp samp )

# @unpack samp = load(datadir("exp_pro", "fulbright123_myega_sampling.jld2"))
# plot(samp, yscale=:log10)


# @time bigsamp = sample(case123_model, NUTS(), MCMCThreads(), 2000, nthreads, initial_params=init_set)
