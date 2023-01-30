using DrWatson;
@quickactivate("IntermittentFields");
include(srcdir("IntermittentFields.jl"));
using .IntermittentFields

Np = 10
dt = 1/(10^2)
N = 2^12
α = 0.25
ξ = 0.6
tmax = 10.
beta = 1.

t, p = doexittime3(Np, dt, N, α, ξ, tmax, beta);

abs.(p[1, :] - p[2, :]) + abs.(p[2, :] - p[3, :]) .> 1
