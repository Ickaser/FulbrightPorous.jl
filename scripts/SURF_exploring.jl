suc_cases = ["2024-06-19-15_MFD_AH.csv",
             "2024-06-21-16_MFD_AH.csv",
             "2024-06-26-11_MFD_AH.csv",
             "2024-06-28-15_MFD_AH.csv",]
man_cases = ["2024-06-06-11_MFD_AH.csv",
             "2024-05-28-12_MFD_AH.csv",
             "2024-05-31-12_MFD_AH.csv",
             "2024-06-04-10_MFD_AH.csv",]

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
for (i,dat) in enumerate(suc_dats)
    # Examine controlled nucleation & whether it worked
    pl1 = blankplot_hrC(xlim=(0, 4), ylim=(-7, 0))
    @df dat exptfplot!(:t, :Tsh_i, label="Tsh", color=:black)
    @df dat exptfplot!(:t, :T1, :T2, :T3, label="")
    # Primary drying temperatures
    pl2 = blankplot_hrC(xlim=(4, 40), ylim=(-40, 0))
    @df dat exptfplot!(:t, :Tsh_i, label="Tsh", color=:black)
    @df dat exptfplot!(:t, :T1, :T2, :T3, label="")
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
    pl1 = blankplot_hrC(xlim=(0, 4), ylim=(-20, 0))
    @df dat exptfplot!(:t, :Tsh_i, label="Tsh", color=:black)
    @df dat exptfplot!(:t, :T1, :T2, :T3, label="")
    # Primary drying temperatures
    pl2 = blankplot_hrC(xlim=(0, 40), ylim=(-40, 0))
    @df dat exptfplot!(:t, :Tsh_i, label="Tsh", color=:black)
    @df dat exptfplot!(:t, :T1, :T2, :T3, label="")
    # Pirani vs CM
    pl3 = @df dat plot(:t, :pch_cm, label="CM", xlim=(0,40), ylim=(0, 200))
    @df dat plot!(:t, :pch_pir, label="Pirani")
    pl = plot(pl1, pl2, pl3, link=:x, layout=(3,1))
    push!(pls, pl)
end
plot(pls..., size=(800, 800))
end
