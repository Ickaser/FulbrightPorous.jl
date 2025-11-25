hf0 = 1u"cm"
l = range(0u"cm", hf0, length=101)
uf = (;unitformat=(x,y)->"$x [$(latexify(y))]")

Rp_cry = RpFormFit(1.0u"cm^2*Torr*hr/g", 16.0u"cm*Torr*hr/g", 0.0u"1/cm")
Rp_amo = RpFormFit(1.0u"cm^2*Torr*hr/g", 16.0u"cm*Torr*hr/g", 1.0u"1/cm")

pl1 = blankplot_lRp(yunit=u"m/s", ylim=(0, Inf), xlim=(0, 1))
plot!(l, Rp_cry.(l), label="crystalline")
plot!(l, Rp_amo.(l), label="amorphous")
plot!(size=(300,400))

pl2 = plot(u"cm", u"s^-1", xlabel=L"h_\mathrm{d}", ylabel=L"\frac{dR_p}{dh_\mathrm{d}}"; uf... )
plot!(l, Rp_cry.A1*(1 .+ Rp_cry.A2*l).^-2, label="crystalline")
plot!(l, Rp_amo.A1*(1 .+ Rp_amo.A2*l).^-2, label="amorphous")
plot!(ylim=(0u"s^-1", Rp_cry.A1*1.1))

plot(pl1, pl2, layout=(2,1), size=(300, 600), left_margin=40Plots.px)

savefig(plotsdir("cryst_v_amorph.svg"))

# ---------------------------
# Typical cycle data

using Latexify

begin
dat1raw = CSV.read(datadir("exp_raw", "Fulbright2025", "2025-01-20-10_MFD_ISW.csv"), Table; header=8)
dat1full = map(x->microfd_fulbright_columnrename(x, dat1raw.CycleTime[1]), dat1raw)
overnightfix(dat1full)
dat1 = filter(x->x.phase==4, dat1full)
dat1.t .-= dat1.t[1]

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

# @df dat1 exptfplot(:t, :T1, :T2, :T3, :T4)

fitdat1 = @df dat1 PrimaryDryFit(:t, 
    (:T1[:t .< 5u"hr"],
    :T2[:t .< 5u"hr"],
    :T3[:t .< 5u"hr"],
    :T4[:t .< 5.9u"hr"],
    );
    t_end=t_end1
)
# plot(fitdat1)

dat2raw = CSV.read(datadir("exp_raw", "Fulbright2025", "2025-01-28-10_MFD_ISW.csv"), Table; header=8)
dat2full = map(x->microfd_fulbright_columnrename(x, dat2raw.CycleTime[1]), dat2raw)
overnightfix(dat2full)
dat2 = filter(x->x.phase==4, dat2full)

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

fitdat2 = @df dat2 PrimaryDryFit(:t, (
    :T1[:t .< 7u"hr"],
    :T2[:t .< 7u"hr"],
    :T3[:t .< 4u"hr"],
    :T4[:t .< 3u"hr"],
    );
    t_end=t_end2
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

fitdat3 = @df dat3 PrimaryDryFit(:t, (
    :T1[:t .< 3u"hr"],
    :T2[:t .< 4u"hr"],
    :T3[:t .< 3.5u"hr"],
    :T4[:t .< 3.5u"hr"],
    );
    t_end=t_end3
)
# plot(fitdat3)
end

begin
Rp2 = RpFormFit(0.002545624558821383u"hr*cm^2*Torr*g^-1", 34.82592716330301u"hr*cm*Torr*g^-1", 13.790538676958578u"cm^-1")
A_p, A_v = π.*get_vial_radii("10R") .^2
fillvol = 2u"mL"
hf0 = fillvol/A_p |> u"cm"
csolid = 0.05u"g/mL"
ρsolution = 1u"g/mL"
pop2 = ParamObjPikal([
    (Rp2, hf0, csolid, ρsolution),
    (Kv, A_v, A_p),
    (pch2, Tsh2),
])
model2 = solve(ODEProblem(pop2), Rodas3())
model2.t .+= ustrip(u"hr", dat2full.t[findfirst(x->x==4, dat2full.phase)])
end

# Latexify.set_default(labelformat=:square)
default(:fontfamily, "sans-serif")

begin
pl1 = @df dat1 exptfplot(:t, :T1, :T2, :T3, :T4, ylabel="Temperature", lw=2)
plot!(Tsh1, tmax=21u"hr", label=L"T_\mathrm{sh}", c=:black, lw=2)
plot!(xlabel="Time", unitformat=:square)
plot!(legend=(0.2, 0.3))
pl2 = twinx()
plot!(pl2, ylabel="Pressure", axiscolor=:green, bordercolor=:green, unitformat=:square)
plot!(pl2, ylim=(0, 180), legend=(0.85, 0.45))
@df dat1[begin+3:end] plot!(pl2, :t, :pch_pir, label="Pirani", c=3, marker=:o, seriestype=:samplemarkers, step=100)
@df dat1[begin+3:end] plot!(pl2, :t, :pch_cm, label="CM", c=3, marker=:square, seriestype=:samplemarkers, step=100)
tendplot!(t_end1, subplot=2)
end

begin
twinx(plot(xwiden=false, size=(450, 280)))
cycledataplot!(dat2full, (:T1, :T2, :T3, :T4), (:Tsh_i), (:pch_pir, :pch_cm), pcolor=:green)
modconvtplot!(model2, lw=3, alpha=0.8)
plot!(subplot=1, legend=(0.65, 0.38), ylim=(-52, 25))
plot!(subplot=2, legend=(0.63, 0.85), ylim=(1, 200))
end
savefig(plotsdir("example_cycle_data.svg"))

l = range(0u"cm", hf0, length=101)
hRps = [calc_hRp_T(pop2, fitdat2; i) for i in 1:4]
markers = [:circle, :square, :diamond, :hexagon]
begin
pl1 = plot(xunit=u"hr", yunit=u"°C", xlabel="t", ylabel="T\\;", unitformat=latexify)
plot!(fitdat2, legend=:bottomright)
pl2 = plot(xunit=u"cm", yunit=u"cm^2*Torr*hr/g", xlabel="h_\\mathrm{d}", ylabel="R_p", unitformat=latexify)
plot!(l, Rp2.(l), label="fit to T", lw=2, c=:black)
for (i, hRp) in enumerate(hRps)
    plot!(hRp[1], hRp[2], label="$i", c=i, marker=markers[i], st=:samplemarkers, step=50, offset=12*(i-1)+1)
end
plot(pl1, pl2, layout=(2,1), size=(340, 320), bottom_margin=-15Plots.px)
end
savefig(plotsdir("example_hRp.svg"))

begin
plot(xticks=[], yticks=[], framestyle=:box, widen=false, legend=false)
plot!(1:3, 1:3, fillrange=0, color=:green, alpha=0.5)
plot!(1:3, [2.5, 1.5, 0.5], fillrange=0, color=:blue, alpha=0.5)
plot!(size=(100, 100))
end
savefig(plotsdir("cartoon_designspace.svg"))