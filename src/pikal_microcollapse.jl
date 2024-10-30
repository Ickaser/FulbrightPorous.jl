const ΔHsub = LyoPronto.ΔHsub
const R = 8.3145u"J/mol/K"
const Mw = 18u"g/mol"

@doc raw"""
"""
function lyo_mc_dae!(du, u, params, t)

    _, hf0, c_solid, ρ_solution = params[1]
    Kshf, Av, Ap, = params[2]
    pch, Tsh = params[3] 
    k0, E, C, S = params[4] 
    
    td = t*u"hr" # Dimensional time
    hf = u[1]*u"cm"
    Tf = u[2]*u"K"
    μ2 = u[3]*u"cm^6"
    if Tf < 0u"K"
        return [NaN, NaN, NaN]
    end

    k = k0*exp(-E/R/T)
    hd = hf0 - hf
    ϕ = (ρ_solution - c_solid)/ρ_solution
    μ1 = hd*ϕ*Ap
    Rp = C * sqrt(R*Tf/Mw) * μ1^2/μ2
    pchl = pch(td)
    Qshf = Av*Kshf(pchl)*(Tsh(td) - Tf)
    Tsub = Tf - Qshf/k_ice/Ap*hf
    dmdt = - Ap*(calc_psub(Tsub)-pch(td))/Rp
    Qsub = uconvert(u"W", dmdt*ΔHsub)

    dhf_dt = min(0u"cm/hr", dmdt/(ρ_solution-c_solid)/Ap) # Cap dhf_dt at 0: no desublimation

    du[1] = ustrip(u"cm/hr", dhf_dt)
    du[2] = ustrip(u"W", Qsub + Qshf)
    du[3] = ustrip(u"cm^6/hr", k*μ1^2 - S(hd)*dhf_dt)
end

const lyo_mc_mm = [1.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 1.0]

@doc raw"""
    lyo_mc_dae_f = ODEFunction(lyo_mc_dae!, mass_matrix=lyo_mc_mm)

Compute the right hand side function for the microcollapse-augmented Pikal model.

The DAE system which is the Pikal model plus pore moment evolution (2 ODEs, one nonlinear algebraic equation for pseudosteady conditions)
is here treated as a constant-mass-matrix implicit ODE system.
The implementation is in [`lyo_mc_dae!`](@ref)

The initial conditions `u0 = [h_f, Tf, μ2]` should be unitless, but are internally assigned to be in `[cm, K, cm^6]`.
The unitless time is taken to be in hours, so derivatives are given in unitless `[cm/hr, K/hr, cm^6/hr]`.

`params` has the form:
```
params = (
    (_, hf0, c_solid, ρ_solution),
    (Kshf, Av, Ap),
    (pch, Tsh) ,
    (k0, E, C, S),
)
```
where those listed following are callables returning `Quantity`s, and the rest are `Quantity`s.
See [`RpFormFit`](@ref LyoPronto.RpFormFit) and [`RampedVariable`](@ref LyoPronto.RampedVariable) for convenience types that can help with the callables.
- `S(hd)` with `hd` a length returns source term for pore size second moment (as a Unitful quantity)
- `Kshf(p)` with `p` a pressure returns heat transfer coefficient (as a Unitful quantity).
- `Tsh(t)`, `pch(t)` return shelf temperature and chamber pressure respectively at time `t`.
"""
lyo_mc_dae_f = ODEFunction(lyo_mc_dae!, mass_matrix=lyo_mc_mm)
