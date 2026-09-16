fp_plot_defaults()
using StatsPlots: groupedbar

begin
Kv = load(datadir("exp_pro", "gravimetric_kv.jld2"))["Kv"]

dat1raw = CSV.read(datadir("exp_raw", "Fulbright2025", "2025-06-06-10_MFD_ISW.csv"), Table; header=8)
dat1nm = map(x->microfd_fulbright_columnrename(x, dat1raw.CycleTime[1]), dat1raw)
dat1 = filter(x->x.phase==4, dat1nm)
dat1.t .-= dat1.t[1]
overnightfix(dat1)

# plot(dat1raw.Phase)

# t_end1 = identify_pd_end(dat1.t, dat1.pch_pir, Val(:der2), window_width=91,  tmin=1u"hr",  tmax=50u"hr")
# t_end1 = identify_pd_end(dat1.t, dat1.pch_pir, Val(:onoff), window_width=91,  tmin=1u"hr",  tmax=50u"hr")
t_end1 = missing

pch1 = RampedVariable(100.0u"mTorr")
@df dat1 plot(:t, :pch_cm, ylim=(0, 200))
@df dat1 plot!(:t, :pch_pir, ylim=(0, 200))
plot!(pch1)
tendplot!(t_end1)

Tsh1 = RampedVariable(uconvert.(u"K", [-45.0, -5]u"°C"), 1.0u"K/minute")
# @df dat1 plot(:t, :Tsh_i, xlim=(0, 8))
# plot!(Tsh1)

pl1 = blankplot_hrC()
@df dat1 exptfplot!(:t, :T1, :T2, :T3, :T4)

fitdat1 = @df dat1 PrimaryDryFit(:t[0u"hr".< :t] , 
    (:T1[0u"hr" .< :t .< 2.9u"hr"],
     :T2[0u"hr" .< :t .< 3.0u"hr"],
     :T3[0u"hr" .< :t .< 1.5u"hr"],
    #  :T4[0u"hr" .< :t .< 3.0u"hr"],
    );
    t_end=t_end1
)
plot(fitdat1)

dat2raw = CSV.read(datadir("exp_raw", "Fulbright2025", "2025-06-09-10_MFD_ISW.csv"), Table; header=8)
dat2nm = map(x->microfd_fulbright_columnrename(x, dat2raw.CycleTime[1]), dat2raw)
dat2 = filter(x->x.phase==4, dat2nm)
dat2.t .-= dat2.t[1]
overnightfix(dat2)

# plot(dat2raw.Phase)

# t_end2 = identify_pd_end(dat2.t, dat2.pch_pir, Val(:2der), window_width=91,  tmin=1u"hr",  tmax=50u"hr")
t_end2 = identify_pd_end(dat2.t, dat2.pch_pir, Val(:onoff), window_width=91,  tmin=1u"hr",  tmax=50u"hr")

pch2 = RampedVariable(100.0u"mTorr")
# @df dat2 plot(:t, :pch_cm, ylim=(0, 200))
# @df dat2 plot!(:t, :pch_pir, ylim=(0, 200))
# plot!(pch2)
# tendplot!(t_end2)

Tsh2 = RampedVariable(uconvert.(u"K", [-45.0,  5]u"°C"), 1u"K/minute")
# @df dat2 plot(:t, :Tsh_i, xlim=(0, 2))
# plot!(Tsh2)

# pl1 = blankplot_hrC()
# @df dat2 exptfplot!(:t, :T1, :T2, :T3, :T4)
# tendplot!(t_end2)

fitdat2 = @df dat2 PrimaryDryFit(:t[0u"hr".< :t] , 
    (:T1[0u"hr" .< :t .< 3u"hr"],
     :T2[0u"hr" .< :t .< 1.9u"hr"],
     :T3[0u"hr" .< :t .< 2.3u"hr"],
    #  :T4[0u"hr" .< :t .< 3u"hr"],
    );
    t_end=t_end2
)
# plot(fitdat2)

dat3raw = CSV.read(datadir("exp_raw", "Fulbright2025", "2025-06-16-10_MFD_ISW.csv"), Table; header=8)
dat3nm = map(x->microfd_fulbright_columnrename(x, dat3raw.CycleTime[1]), dat3raw)
dat3 = filter(x->x.phase==4, dat3nm)
dat3.t .-= dat3.t[1]
overnightfix(dat3nm)
overnightfix(dat3)

# plot(dat3raw.Phase)

# t_end3 = identify_pd_end(dat3.t, dat3.pch_pir, Val(:2der), window_width=91,  tmin=1u"hr",  tmax=50u"hr")
t_end3 = identify_pd_end(dat3.t, dat3.pch_pir, Val(:onoff), window_width=91,  tmin=1u"hr",  tmax=50u"hr")

pch3 = RampedVariable(100.0u"mTorr")
# @df dat3 plot(:t, :pch_cm, ylim=(0, 200))
# @df dat3 plot!(:t, :pch_pir, ylim=(0, 200))
# plot!(pch3)
# tendplot!(t_end3)

Tsh3 = RampedVariable(uconvert.(u"K", [-45.0,  -10]u"°C"), 1u"K/minute")
# @df dat3 plot(:t, :Tsh_i, xlim=(0, 2))
# plot!(Tsh3)

pl1 = blankplot_hrC()
# @df dat3 exptfplot!(:t, :T1, :T2, :T3, :T4)
# tendplot!(t_end3)

fitdat3 = @df dat3 PrimaryDryFit(:t[0u"hr".< :t] , 
    (
    :T1[0u"hr" .< :t .< 12u"hr"],
    :T2[0u"hr" .< :t .< 12u"hr"],
    :T3[0u"hr" .< :t .< 3.5u"hr"],
    :T4[0u"hr" .< :t .< 5u"hr"],
    );
    t_end=t_end3
)
# plot(fitdat3)

dat4raw = CSV.read(datadir("exp_raw", "Fulbright2025", "2025-06-18-10_MFD_ISW.csv"), Table; header=8)
dat4nm = map(x->microfd_fulbright_columnrename(x, dat4raw.CycleTime[1]), dat4raw)
dat4 = filter(x->x.phase==4, dat4nm)
dat4.t .-= dat4.t[1]
overnightfix(dat4nm)
overnightfix(dat4)

# plot(dat4raw.Phase)

# t_end4 = identify_pd_end(dat4.t, dat4.pch_pir, Val(:2der), window_width=91,  tmin=1u"hr",  tmax=50u"hr")
t_end4 = identify_pd_end(dat4.t, dat4.pch_pir, Val(:onoff), window_width=91,  tmin=5u"hr",  tmax=50u"hr")

pch4 = RampedVariable(100.0u"mTorr")
@df dat4 plot(:t, :pch_cm, ylim=(0, 200))
@df dat4 plot!(:t, :pch_pir, ylim=(0, 200))
plot!(pch4)
tendplot!(t_end4)

Tsh4 = RampedVariable(uconvert.(u"K", [-45.0,  0]u"°C"), 1u"K/minute")
# @df dat4 plot(:t, :Tsh_i, xlim=(0, 2))
# plot!(Tsh4)

# pl1 = blankplot_hrC()
# @df dat4 exptfplot!(:t, :T1, :T2, :T3, :T4)
# tendplot!(t_end4)

fitdat4 = @df dat4 PrimaryDryFit(:t[0u"hr".< :t] , 
    (
    :T1[0u"hr" .< :t .< 7u"hr"],
    :T2[0u"hr" .< :t .< 7u"hr"],
    :T3[0u"hr" .< :t .< 7u"hr"],
    :T4[0u"hr" .< :t .< 2u"hr"],
    );
    t_end=t_end4
)
# plot(fitdat4, nmarks=10, showline=true)


end

blankplot_hrC()
plot!(fitdat1, c=1)
plot!(fitdat2, c=2)
plot!(fitdat3, c=3)
plot!(fitdat4, c=4)

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
k0 = 0.0u"1/s"

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
# @reset pop2.Rp.R0 = 1e-8u"cm^2*Torr*hr/g"
# @reset pop2.Rp.A1 = 32.3u"cm*Torr*hr/g"
pomc2 = @set pomc1.Tsh = Tsh2

pop3 = @set pop1.Tsh = Tsh3
@reset pop3.hf0 = 3u"mL"/A_p |> u"cm"
pop4 = @set pop3.Tsh = Tsh4
pomc3 = @set pomc1.Tsh = Tsh3
@reset pomc3.hf0 = 3u"mL"/A_p |> u"cm"
pomc4 = @set pomc3.Tsh = Tsh4
end

@df dat1nm exptfplot(:t, :T1, :T2, :T3, :T4, xlim=(1,2))
@df dat2nm exptfplot(:t, :T1, :T2, :T3, :T4, xlim=(1,2))
@df dat3nm exptfplot(:t, :T1, :T2, :T3, :T4, xlim=(1,2))
@df dat4nm exptfplot(:t, :T1, :T2, :T3, :T4, xlim=(1,2))

# blankplot_hrC()
# modconvtplot!(solve(ODEProblem(pop2), LyoPronto.odealg_chunk2))

cases = Table([ (name="1", pop=pop1, pomc=pomc1, fitdat=fitdat1, dat=dat1),
                (name="2", pop=pop2, pomc=pomc2, fitdat=fitdat2, dat=dat2),
                (name="3", pop=pop3, pomc=pomc3, fitdat=fitdat3, dat=dat3),
                (name="4", pop=pop4, pomc=pomc4, fitdat=fitdat4, dat=dat4),
])

trans_Rp = Rp_transform_basic(1.4u"cm^2*Torr*hr/g", 16.0u"cm*Torr*hr/g", 1e-2u"1/cm")
obj_Rpf = OptimizationFunction((x,y)->obj_pd(x,y,tweight=1e-1), AutoForwardDiff(chunksize=3))
@time obj_Rpf([0.0, 1.0, 0.0], (trans_Rp, pop1, fitdat1))

single_Rps = map(cases) do case
    (;name, pop, fitdat, dat) = case
    u0 = [-4.0, 1.0, 0.0]
    @time opt = solve(OP(obj_Rpf, u0, (trans_Rp, pop, fitdat)), optalg)#, g_tol=1e-6)
    Rp = transform(trans_Rp, opt.u).Rp
    sol = gen_sol_pd(opt.u, trans_Rp, pop, fitdat)
    @info "analyis" name opt Rp sol.retcode
    hRps = map(x->calc_hRp_T(case.pop, case.fitdat, i=x), 1:length(fitdat.Tfs))
    return (opt=opt, sol_p=sol, Rp=Rp, hRps=hRps)
end;

obj_n = OptimizationFunction((x,y)->LyoPronto.objn_pd(x,y,tweight=1e-1), AutoForwardDiff())
rp_all = @time solve(OP(obj_n, fill(0.0, 3), (trans_Rp, cases.pop, cases.fitdat)), optalg;)

cases_Rps = Table(cases, single_Rps)

errs = Table(map(cases_Rps) do case
    shared = gen_sol_pd(rp_all.u, trans_Rp, case.pop)
    PP = case.sol_p
    errf = (model, tweight)->obj_expT(model, case.fitdat; tweight)
    sh_T = sqrt(errf(shared, 0))*u"K"
    sep_T = sqrt(errf(PP, 0))*u"K"
    t_end_est = ismissing(case.fitdat.t_end) ? NaN : ustrip(u"hr", sum(case.fitdat.t_end)/2)
    sh_t = sqrt(errf(shared, 1) - errf(shared, 0))/t_end_est
    sep_t = sqrt(errf(PP, 1) - errf(PP, 0))/t_end_est
    return (;name=case.name, sh_T, sep_T, sh_t, sep_t)
end)

begin
e1 = @df errs groupedbar(:name, [:sh_T :sep_T], label="", ylabel=L"$T_\mathrm{f}$ error", c=[2 4], lw=1)
e2 = @df errs groupedbar(:name, [:sh_t :sep_t], label="", ylabel=L"$t$ rel. error", c=[2 4], lw=1)
plot!(e1, left_margin=30Plots.px)
plot!(e2, left_margin=30Plots.px)
vspan!(e2, [0, 1], c=:lightgray, label="")
pls = map(cases_Rps) do case 
    plT = blankplot_hrC()
    @df case.dat exptfplot!(:t, :T1, :T2, :T3, :T4, label="", seriestype=:line, marker=:none)
    plot!(case.fitdat, title=case.name, nmarks=10, showline=true, markeralpha=0.6)
    modconvtplot!(case.sol_p, label="Pikal fit to case", alpha=0.7)
    modconvtplot!(gen_sol_pd(rp_all.u, trans_Rp, case.pop), c=:purple, label="Pikal fit to all", alpha=0.7)
    plot!(ylim=(-45, 10))
    plot!(legend=:none)
    if case.name == "1"
        plot!(legend=(-0.2, -0.6), legendfontsize=14, left_margin=20Plots.px)
    else
        plot!(ylabel="", left_margin=-15Plots.px)
    end
    plot!(case.pop.Tsh, tmax=case.fitdat.t[end], c=:black, label=L"T_\mathrm{sh}")
    return plT
end
plot(pls..., e1, e2, layout=@layout([grid(1,4); _ grid(1,2){0.8w}]), size=(800, 400))
end
savefig(plotsdir("fulbright_dextran_T.svg"))
# Manually adjust legend before converting to PDF

l = range(0u"cm", hf0, length=100)
begin
marks = [:o, :utriangle, :square, :diamond]
pls = map(cases_Rps) do case
    plr = blankplot_lRp(title=case.name)
    for (i, hRp) in enumerate(case.hRps)
        hdT, RpT = hRp
        plot!(hdT, RpT, c=:gray, marker=marks[i], seriestype=:samplemarkers, step=40, offset=mod1(i, 4)*10, alpha=0.8, label="Rp from T, $i")
    end
    plot!(l, case.Rp.(l), c=1, label="3PP fit directly to T")
    plot!(l, transform(trans_Rp, rp_all.u).Rp.(l), c=2, label="3PP fit to all T")
    plot!(legend = :topleft)
    return plr
end
plot(pls..., layout=(2,2), size=(800, 600))
end

begin
marks = [:circle, :utriangle, :square, :diamond]
csamp = [90, 80, 70, 50]
cs = stack([colormap("blues")[csamp],
            colormap("oranges")[csamp],
            colormap("greens")[csamp],
            palette([:purple, :white], 7)[1:4],
            ])

pls = map((dat1nm, dat2nm, dat3nm, dat4nm), 1:4) do dat, i
    pln = blankplot_hrC(ylabel=L"T")
    labs = ["Exp. $i, \$T_\\mathrm{f$j}\$" for j in 1:4] |> permutedims
    nucdat = filter(x->x.phase==1, dat)
    if i ∈ 3:4
        @df nucdat exptfplot!(:t, :T1, :T2, :T3, :T4, labels=labs, showline=true, nmarks=20, c=cs[:,i]', markershape=permutedims(marks))
    else
        @df nucdat exptfplot!(:t, :T1, :T2, :T3, labels=labs, showline=true, nmarks=20, c=cs[:,i]', markershape=permutedims(marks))
    end
    i ∈ 1:3 && plot!(xlabel="", bottom_margin=-15Plots.px)
    i ∈ 1:3 ? plot!(xticks=(1.0:0.5:2.0, [])) : plot!(xticks=(1.0:0.5:2.5))
    plot!(xlim=(1,2.06), ylim=(-17, 1), legend=(1.2, 0.8))
end
pln = plot(pls..., layout=(4,1), size=(500,800))

plr = blankplot_lRp()
for case in cases_Rps[1:4]
    j = tryparse(Int, case.name)
    for (i, hRp) in enumerate(case.hRps)
        hdT, RpT = hRp
        plot!(hdT, RpT, c=cs[i,j], marker=marks[i], 
        seriestype=:samplemarkers, step=60, offset=mod1(i, 4)*8+ (mod(j,0:3))*10, 
        linealpha=0.3, 
        # label=(i==1 ? "Rp from T, $(case.name)" : nothing)
        label=""
        )
    end
end
plot!(xlim=(0, 0.3), ylim=(0, 5))

plot!(pln, plr, layout=@layout([a _ b]), size=(600, 500))
end
savefig(plotsdir("fulbright_dextran_rawRp.svg"))
savefig(plotsdir("fulbright_dextran_rawRp.pdf"))

function regress_line(xs, ys)
    n = length(xs)
    @assert n == length(ys)
    sx = sum(xs)
    sy = sum(ys)
    sxx = sum(xs.^2)
    syy = sum(ys.^2)
    sxy = sum(xs .* ys)
    denom = n*sxx - sx^2
    m = (n*sxy - sx*sy) / denom
    b = (sxx*sy - sx*sxy) / denom
    return m, b
end


begin
ends = fill(NaN*u"cm", (4, 4))
ends[1,1] = 0.16u"cm"
ends[1,3] = 0.06u"cm"
ends[2,2] = 0.12u"cm"
ends[4,1] = 0.19u"cm"
ends[4,2] = 0.21u"cm"
ends[4,3] = 0.21u"cm"
ends[3,1] = 0.19u"cm"
ends[3,2] = 0.25u"cm"
ends[3,3] = 0.16u"cm"

lines = map(enumerate(cases_Rps[1:4])) do (j, case)
    Table(map(enumerate(case.hRps)) do (i, hRp)
        hdT, RpT = hRp
        if !isnan(ends[j,i])
            hdT = filter(x->x<ends[j,i], hdT)
            RpT = RpT[begin:length(hdT)]
        end
        A1, R0 = regress_line(hdT, RpT) 
        (;R0, A1)
    end)
end
pls = map(cases_Rps[1:4]) do case
    plr = blankplot_lRp()
    j = tryparse(Int, case.name)
    for (i, hRp) in enumerate(case.hRps)
        hdT, RpT = hRp
        if !isnan(ends[j,i])
            hdT = filter(x->x<ends[j,i], hdT)
            RpT = RpT[begin:length(hdT)]
        end
        (;R0, A1) = lines[j][i]
        plot!(l, l*A1.+R0, c=cs[i,j], label="")

        plot!(hdT, RpT, c=cs[i,j], marker=marks[i], 
        seriestype=:samplemarkers, step=10, offset=mod1(i, 4)*8+ (mod(j,0:3))*10, 
        linealpha=0.0, 
        label=(i==1 ? "Rp from T, $(case.name)" : nothing)
        # label=""
        )
    end
    plot!(xlim=(0.0, 0.2), ylim=(0, 3))
end
plot(pls...)


end

begin
plot(layout=(4,1), size=(400, 400))
t_res = map((dat1nm, dat2nm, dat3nm, dat4nm), 1:4) do dat, i
    labs = ["Exp. $i, \$T_\\mathrm{f$j}\$" for j in 1:4] |> permutedims
    Ts = map(i->getproperty(dat, Symbol("T$i"))[begin:180], 1:4)
    i_nucs = map(Ts) do T
        i_nuc = argmax(savitzky_golay(T.|>u"K", 3, 1; deriv=1).y) -1
        i_nuc
    end
    i_ress = map(Ts, i_nucs) do T, i_nuc
        findlast(x->x > -10u"degC", T)
    end

    nucdat = filter(x->x.t < 3u"hr", dat)
    @df nucdat exptfplot!(:t, :T1, :T2, :T3, :T4, labels=labs, showline=true, nmarks=20, c=cs[:,i]', markershape=permutedims(marks), subplot=i)
    for (j, (i_nuc, i_res)) in enumerate(zip(i_nucs, i_ress))
        # plot!([dat.t[i_nuc], dat.t[i_res]], Ts[j][[i_nuc,i_res]], c=:black, arrow=true, label="", subplot=i)
        plot!(dat.t[i_nuc:i_nuc], Ts[j][i_nuc:i_nuc], marker=:star, c=:red, label="", subplot=i)
    end
    i ∈ 1:3 && plot!(xlabel="", bottom_margin=-15Plots.px)
    plot!(xlim=(1,2.8), ylim=(-15, 0))
    return (i_ress .- i_nucs).*u"minute"
end
tres_plot = plot!()
end

begin
pl_a1t = plot(u"hr", u"hr*cm*Torr/g", size=(400, 300))
for (tt, line, i) in zip(t_res, lines, 1:4)
    scatter!(tt[1:length(line)], line.A1, c=cs[:,i], marker=marks, label="")
end
plot!(xlabel="Time after nucleation \n above -10°C", ylabel="Initial slope", unitformat=:square)
plot!(xticks=0.2:0.1:0.6, xlim=(0.2, 0.54))
end

begin
plot(pln, plr, pl_a1t, layout=@layout([a b c]), size=(1000, 400))
plot!(left_margin=20Plots.px, subplot=1)
for i in 4:6; plot!(subplot=i, bottom_margin=30Plots.px); end
for i in 1:3; plot!(subplot=i, xticks=(1.0:0.5:2.0, []), right_margin=150Plots.px); end
plot!(subplot=5, right_margin=20Plots.px, xticks=0:0.1:0.3)
plot!()
end
savefig(plotsdir("fulbright_dextran_A1_vs_crysttime.svg"))
# Fix the legend before exporting to PDF
