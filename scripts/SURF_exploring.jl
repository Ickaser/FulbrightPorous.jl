suc_cases = ["2024-06-19-15_MFD_AH.csv",
             "2024-06-21-16_MFD_AH.csv",
             "2024-06-26-11_MFD_AH.csv",
             "2024-06-28-15_MFD_AH.csv",]
man_cases = ["2024-06-06-11_MFD_AH.csv",
             "2024-05-28-12_MFD_AH.csv",
             "2024-05-31-12_MFD_AH.csv",
             "2024-06-04-10_MFD_AH.csv",]

function microfd_columnrename(row, t0 = nothing)
    # nt = (tstamp = DateTime(row.Timestamp, dateformat"mm/dd/yyyy H:M:S"),
    if !isnothing(t0)
        t = uconvert(u"hr", row.CycleTime - t0)
        if t < 0u"hr"
            t += 24u"hr"
        end
    else
        t = nothing
    end
    nt = (tstamp = row.CycleTime,
          t = t,
          phase = row.Phase,
          step = row.Step,
          pch_sp = row.VacSetPT*u"mTorr",
          pch_pir = row.VacPirani*u"mTorr",
          pch_cm = row.VacCPM*u"mTorr",
          Tsh_sp = row.ShelfSetPT*u"°C",
          Tsh_i = row.ShelfInlet*u"°C",
          T1 = row.TP1*u"°C",
          T2 = row.TP2*u"°C",
          T3 = row.TP4*u"°C", # Not a mistake: quirk of experimental apparatus that connection TP3 was more annoying than TP4
          # if a column has names with whitespace 
          # or special characters in it, do this:
          # mfr = row.var"Mass Flow Rate"
    )
    return nt
end
function readdat(case)
    tab = CSV.read(datadir("exp_raw", "SURF2024", case), Table, header=8)
    tab = map(row -> microfd_columnrename(row, tab[1].CycleTime), tab)
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
for dat in suc_dats
    # Examine controlled nucleation & whether it worked
    pl1 = blankplot_hrC(xlim=(0, 4), ylim=(-7, 0))
    plot!(dat.t, dat.Tsh_i, label="")
    plot!(dat.t, dat.T1, label="")
    plot!(dat.t, dat.T2, label="")
    plot!(dat.t, dat.T3, label="")
    # Primary drying temperatures
    pl2 = blankplot_hrC(xlim=(4, 40), ylim=(-40, 0))
    plot!(dat.t, dat.Tsh_i, label="")
    plot!(dat.t, dat.T1, label="")
    plot!(dat.t, dat.T2, label="")
    plot!(dat.t, dat.T3, label="")
    hline!([-32], label="")
    # Pirani vs CM
    pl3 = plot(dat.t, dat.pch_cm,label="CM", xlim=(0,40), ylim=(0, 200))
    plot!(dat.t, dat.pch_pir, label="Pirani")
    pl = plot(pl1, pl2, pl3, link=:x, layout=(3,1))
    push!(pls, pl)
end
plot(pls..., size=(800, 800))
end

begin
pls = []
for dat in man_dats
    # Examine controlled nucleation & whether it worked
    pl1 = blankplot_hrC(xlim=(0, 4), ylim=(-20, 0))
    plot!(dat.t, dat.Tsh_i, label="")
    plot!(dat.t, dat.T1, label="")
    plot!(dat.t, dat.T2, label="")
    plot!(dat.t, dat.T3, label="")
    # Primary drying temperatures
    pl2 = blankplot_hrC(xlim=(0, 40), ylim=(-40, 0))
    plot!(dat.t, dat.Tsh_i, label="")
    plot!(dat.t, dat.T1, label="")
    plot!(dat.t, dat.T2, label="")
    plot!(dat.t, dat.T3, label="")
    # Pirani vs CM
    pl3 = plot(dat.t, dat.pch_cm,label="CM", xlim=(0,40), ylim=(0, 200))
    plot!(dat.t, dat.pch_pir, label="Pirani")
    pl = plot(pl1, pl2, pl3, link=:x, layout=(3,1))
    push!(pls, pl)
end
plot(pls..., size=(800, 800))
end
