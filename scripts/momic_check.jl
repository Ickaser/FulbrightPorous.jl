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


