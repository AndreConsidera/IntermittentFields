using DrWatson;
@quickactivate("Intermittency Paradox")
include(srcdir("intermittentfields_mod.jl"))
using .IntermittentFields
using Plots

N = 2^15;
eta = 1/N;
ξ = 0.8
α = 0.2
g1 = UnitaryWhiteNoise(N, true);
g2 = UnitaryWhiteNoise(N, true);

# Gaussian
u, r = quickrealization(N, eta, ξ);
plot(r, u)

u0, r = quickrealization(N, eta, ξ, g1);
plot(r, u0)
# intermittent
u, r = quickrealization(N, eta, ξ, α);
plot(r, u)

u, r = quickrealization(N, eta, ξ, α, g1);
plot(r, u)

u1, r = quickrealization(N, eta, ξ, α, g1, g2);
plot(r, u1)

plot(r,[u0,u1], ylims=(-6,6))



