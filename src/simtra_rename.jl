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