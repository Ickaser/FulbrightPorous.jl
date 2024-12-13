using QuadGK
using Distributions

# Choices sampled from my head and GitHub Copilot
all_distributions = [
    Beta(1, 3),
    Beta(2, 2),
    Normal(0, 1),
    Normal(10, 5),
    Exponential(1),
    Exponential(10),
    Gamma(1, 1),
    Gamma(2, 2),
    Uniform(-1, 1),
    Uniform(-10, 10),
    Uniform(0, 1),
]

# For each distribution, integrate the PDF with a range of moment orders

function moment_integral(dist, n)
    f(x) = x^n * pdf(dist, x)
    _, b = extrema(dist)
    a = 0
    integral, _ = quadgk(f, a, b)
    return integral
end

tab = Table(
    map(all_distributions) do dist 
    # println("Distribution: ", dist)
    name = replace(string(dist), r"\{Float64\}"=>"")
    ords = [0, 1, 2, 1.5]
    moms = [moment_integral(dist, n) for n in ords]
    mom15 = moms[1]^(-1//8)*moms[2]^(3//4)*moms[3]^(3//8)
    return (dist=name, m0=moms[1], m1=moms[2], m2=moms[3], m32ex=moms[4], m32est=mom15)
end
)

tab

plot([1e-1, 1e2], [1e-1,1e2], color=:black, lw=2, label="parity")
for row in tab
    scatter!([row.m32ex], [row.m32est], label=row.dist)
end
plot!(xlabel="Exact 1.5th moment", ylabel="Interpolated 1.5th moment", legend=:topleft, scale=:log10)
savefig(plotsdir("momic_check.svg"))
savefig(plotsdir("momic_check.pdf"))