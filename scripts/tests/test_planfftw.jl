using DrWatson
@quickactivate("IntermittentFields");
include(srcdir("IntermittentFields.jl"))
using .IntermittentFields
using FFTW
using Plots;plotly()

# constructing kernel
N=2^15
r = Array(range(-pi,stop=pi,length=N),); 
eta = 0.0;
ξ = 2/3
expker = CovarianceKernel(r, eta, ξ);
logker = SingularKernel(r, eta);

# testing different methods of realization function
g1 = UnitaryWhiteNoise(div(N, 2) + 1);
Pinv = plan_irfft(copy(g1.wk), length(expker.r); flags = FFTW.MEASURE, timelimit = Inf);
u = realization(expker, g1, Pinv)
v = realization(expker, g1)
plot(r, [u, v])

g1 = UnitaryWhiteNoise(div(N, 2) + 1);
g2 = UnitaryWhiteNoise(div(N, 2) + 1);
α = 0.25
Γ = GmcNoise(logker,g1,g2,α);
u = realization(expker, Γ);
v = realization(expker, Γ, Pinv)
plot(r, [u, v])

# testing GmcNoise
Pforward = plan_rfft(copy(expker.r); flags = FFTW.MEASURE, timelimit = Inf);
Γ = GmcNoise(logker, g1, g2, α, Pforward, Pinv);
γ = GmcNoise(logker, g1, g2, α);
plot(r, [Γ.wr, γ.wr])

# testing exit times
Np = 200
dt = 1/(2*10^2)
N = 2^12
α = 0.25
ξ = 0.6
tmax = 10.
beta = 1.
planning_effort = FFTW.MEASURE

using Random;
Random.seed!(12345678);
et1 = doexittime(Np, dt, N, α, ξ, tmax, beta, planning_effort);
Random.seed!(12345678);
et2 = doexittime(Np, dt, N, α, ξ, tmax, beta);
et1 == et2

# time
Random.seed!(12345678);
@time doexittime(Np, dt, N, α, ξ, tmax, beta, planning_effort);
Random.seed!(12345678);
@time doexittime(Np, dt, N, α, ξ, tmax, beta);