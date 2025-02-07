# Double check this math against the data from MicroFD, might be off
Kv = load(datadir("exp_pro", "gravimetric_kv.jld2"))["Kv"]

dat1raw = CSV.read(datadir("exp_raw", "Fulbright2025", "2025-01-20-10_MFD_ISW.csv"), Table; header=8)
dat1pd = filter(x->x.Phase==4, dat1raw)
dat1 = map(x->microfd_fulbright_columnrename(x, dat1pd.CycleTime[1]), dat1pd)
for i in eachindex(dat1)[2:end]
    if dat1.t[i] < dat1.t[i-1]
        dat1.t[i:end] .+= 24u"hr"
    end
end
plot(dat1raw.Phase)

pch1_pir_der2 = savitzky_golay(dat1.pch_pir, 91, 3, deriv=2, rate=1/u"minute").y
t_end1 = dat1.t[argmax(pch1_pir_der2[1u"hr" .< dat1.t .< 50u"hr" ])] + 1u"hr"

pch1 = RampedVariable(100u"mTorr")
@df dat1 plot(:t, :pch_cm, ylim=(0, 200))
@df dat1 plot!(:t, :pch_pir, ylim=(0, 200))
plot!(pch1)
tendplot!(t_end1)

Tsh1 = RampedVariable(uconvert.(u"K", [-47, -20]u"°C"), 0.5u"K/minute")
@df dat1 plot(:t, :Tsh_i, xlim=(0, 2))
plot!(Tsh1)

@df dat1 exptfplot(:t, :T1, :T2, :T3, :T4)

fitdat1 = @df dat1 PrimaryDryFit(:t, (
    :T1[:t .< 5u"hr"],
    :T2[:t .< 5u"hr"],
    :T3[:t .< 5u"hr"],
    :T4[:t .< 5.9u"hr"],
    ),
    t_end1
)
plot(fitdat1)

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
@df dat2 plot(:t, :pch_cm, ylim=(0, 200))
@df dat2 plot!(:t, :pch_pir, ylim=(0, 200))
plot!(pch2)
tendplot!(t_end2)

Tsh2 = RampedVariable(uconvert.(u"K", [-47, -10]u"°C"), 0.5u"K/minute")
@df dat2 plot(:t, :Tsh_i, xlim=(0, 2))
plot!(Tsh2)

@df dat2 exptfplot(:t, :T1, :T2, :T3, :T4)

fitdat2 = @df dat2 PrimaryDryFit(:t, (
    :T1[:t .< 7u"hr"],
    :T2[:t .< 7u"hr"],
    :T3[:t .< 5u"hr"],
    :T4[:t .< 3u"hr"],
    ),
    t_end2
)
plot(fitdat2)

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
@df dat3 plot(:t, :pch_cm, ylim=(0, 200))
@df dat3 plot!(:t, :pch_pir, ylim=(0, 200))
plot!(pch3)
tendplot!(t_end3)

Tsh3 = RampedVariable(uconvert.(u"K", [-47, 0]u"°C"), 0.5u"K/minute")
@df dat3 plot(:t, :Tsh_i, xlim=(0, 2))
plot!(Tsh3)

@df dat3 exptfplot(:t, :T1, :T2, :T3, :T4)

fitdat3 = @df dat3 PrimaryDryFit(:t, (
    :T1[:t .< 3u"hr"],
    :T2[:t .< 4u"hr"],
    :T3[:t .< 3.5u"hr"],
    :T4[:t .< 3.5u"hr"],
    ),
    t_end3
)
plot(fitdat3)

blankplot_hrC()
plot!(fitdat1, c=1)
plot!(fitdat2, c=2)
plot!(fitdat3, c=3)




# ---------------
# blankplot_hrC(xlim=(0, 10), ylim=(-40, -30))
# @df dat1 exptfplot!(:t, :T1, :T2, :T3, :T4, c=1, labsuffix=", -20°C")
# @df dat2 exptfplot!(:t, :T1, :T2, :T3, :T4, c=2, labsuffix=", -10°C")
# savefig(plotsdir("sucrose_2cases.png"))