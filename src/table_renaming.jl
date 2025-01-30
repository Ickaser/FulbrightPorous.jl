export simtra_rename
 
repmissing(x) = ismissing(x) ? NaN : x
function simtra_rename(row)

    nt = (t = row.var"Delta Time (h)"*u"hr",
        Tsh_i = row.var"SHELF IN"*u"°C",
        Tsh_o = row.var"SHELF OUT"*u"°C",
        Tsh = row.var"SHELF SETPT"*u"°C",
        pch_pir = row.PIRANI*u"mTorr",
        pch_cm = row.CM*u"mTorr",
        T1 = repmissing(row.TP01)*u"°C",
        T2 = repmissing(row.TP02)*u"°C",
        T3 = repmissing(row.TP03)*u"°C",
        mfr = repmissing(row.var"MFR (g/vial/hr)")*u"g/hr",
        Rp = repmissing(row.Rp)*u"cm^2*hr*Torr/g",
    )
    return nt
end

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