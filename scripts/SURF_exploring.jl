suc_cases = ["2024-06-19-15_MFD_AH.csv",
             "2024-06-21-16_MFD_AH.csv",
             "2024-06-26-11_MFD_AH.csv",
             "2024-06-28-15_MFD_AH.csv",]
man_cases = ["2024-06-06-11_MFD_AH.csv",
             "2024-05-28-12_MFD_AH.csv",
             "2024-05-31-12_MFD_AH.csv",
             "2024-06-04-10_MFD_AH.csv",]
fp_plot_defaults()
nmarks = 15

function readdat(case)
    tab = CSV.read(datadir("exp_raw", "SURF2024", case), Table, header=8)
    tab = map(row -> microfd_SURF_columnrename(row, tab[1].CycleTime), tab)
    for i in 2:length(tab)
        if tab.t[i] < tab.t[i-1]
            tab.t[i:end] .+= 24u"hr"
        end
    end
    tab
end

suc_dats = map(readdat, suc_cases);
man_dats = map(readdat, man_cases);

# Examine process data in wide view
begin
pls = []
for (i,dat) in enumerate(suc_dats)
    # Examine controlled nucleation & whether it worked
    pl1 = blankplot_hrC(xlim=(0, 4), ylim=(-7, 0))
    @df dat exptfplot!(:t, :Tsh_i, label="Tsh", color=:black, nmarks=1, sampmarks=true)
    @df dat exptfplot!(:t, :T1, :T2, :T3, label="", nmarks=300)
    # Primary drying temperatures
    pl2 = blankplot_hrC(xlim=(4, 40), ylim=(-40, 0))
    @df dat exptfplot!(:t, :Tsh_i, label="Tsh", color=:black, nmarks=1, sampmarks=true)
    @df dat exptfplot!(:t, :T1, :T2, :T3, label="", nmarks=30)
    hline!([-32], label="")
    # Pirani vs CM
    pl3 = @df dat plot(:t, :pch_cm, label="CM", xlim=(0,40), ylim=(0, 200))
    @df dat plot!(:t, :pch_pir, label="Pirani")
    # Assemble
    pl = plot(pl1, pl2, pl3, link=:x, layout=(3,1))
    push!(pls, pl)
end
plot(pls..., size=(800, 800))
end

begin
pls = []
for (i,dat) in enumerate(man_dats)
    # Examine controlled nucleation & whether it worked
    pl1 = blankplot_hrC(xlim=(0, 4), ylim=(-20, 3))
    @df dat exptfplot!(:t, :Tsh_i, label="Tsh", color=:black, nmarks=1, sampmarks=true)
    @df dat exptfplot!(:t, :T1, :T2, :T3, label=""; nmarks=300)
    # Primary drying temperatures
    pl2 = blankplot_hrC(xlim=(0, 40), ylim=(-40, 0))
    @df dat exptfplot!(:t, :Tsh_i, label="Tsh", color=:black; nmarks=1, sampmarks=true)
    @df dat exptfplot!(:t, :T1, :T2, :T3, label=""; nmarks=30)
    # Pirani vs CM
    pl3 = @df dat plot(:t, :pch_cm, label="CM", xlim=(0,40), ylim=(0, 200))
    @df dat plot!(:t, :pch_pir, label="Pirani")
    pl = plot(pl1, pl2, pl3, link=:x, layout=(3,1))
    push!(pls, pl)
end
plot(pls..., size=(800, 800))
end


Kv_transform = LyoPronto.K_transform_basic(20.0u"W/m^2/K")
Rp_transform = LyoPronto.Rp_transform_basic(1.0u"cm^2*Torr*hr/g", 16.0u"cm*Torr*hr/g", 0.1u"cm^-1")

transform_4R1K = as((;separate=as(Vector,  Rp_transform, 4),
              shared=Kv_transform,))

man_times = [[14u"hr", 15u"hr", 14u"hr"], 
             [5u"hr", 16u"hr", 15u"hr"], 
             [10u"hr", 10u"hr", 13u"hr"], 
             [14u"hr", 14u"hr", 15u"hr"],]
man_pdfit = map(man_dats, man_times) do dat, times
    dat_pd = filter(row -> row.phase == 4, dat)
    dat_pd.t .-= dat_pd.t[1]
    t_end = identify_pd_end(dat_pd.t, dat_pd.pch_pir, Val(:onoff); tmin=2u"hr")
    @df dat_pd PrimaryDryFit(:t, (:T1[:t.<times[1]], :T2[:t.<times[2]], :T3[:t.<times[3]]); t_end)
end;
plot(map(x->plot(x, nmarks=15, legend=:bottomright), man_pdfit)..., )


begin
rad_i, rad_o = get_vial_radii("10R")
Ap = π*rad_i^2
Av = π*rad_o^2
fillvol = 3u"mL"
hf0 = fillvol/Ap

pch = RampedVariable(100u"mTorr")
Tsh = RampedVariable([-40, -10]u"°C" .|>u"K", 1u"K/minute")

Rp_g = RpFormFit(0.8u"cm^2*Torr*hr/g", 14u"cm*Torr*hr/g", 1u"1/cm")
Kshf_g = ConstPhysProp(15.0u"W/m^2/K")

csolid = 0.05u"g/mL"
ρsolution = 1u"g/mL"

pop = ParamObjPikal([
    (Rp_g, hf0, csolid, ρsolution),
    (Kshf_g, Av, Ap),
    (pch, Tsh),
])
end

man_pdfit


obj_npf = OptimizationFunction((x,y)->objn_pd(x,y,tweight=1e-2), AutoForwardDiff(chunksize=5))
obj_npf(fill(0.0, 13), (transform_4R1K, pop, man_pdfit))

OP = OptimizationProblem(obj_npf, fill(0.0, 13), (transform_4R1K, pop, man_pdfit))
opt = solve(OP, optalg)

params = transform(transform_4R1K, opt.u)
sols = gen_nsol_pd(opt.u, transform_4R1K, pop, man_pdfit);

Kshf = params.shared.Kshf.val |> u"cal/s/cm^2/K"
Rps = Table(Table(params.separate).Rp)
@info "Fitted parameters:" Kshf Rps

plot(map(sols, man_pdfit) do sol, fit
    blankplot_hrC()
    plot!(fit; nmarks=15, )
    modconvtplot!(sol; label="model")
end...)