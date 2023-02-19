using DrWatson
@quickactivate("IntermittentFields");
using FFTW
using Distributions
using Random
using Statistics
using Interpolations
include(srcdir("norms.jl"));
include(srcdir("correlations.jl"));
include(srcdir("kernels.jl"));
include(srcdir("noise.jl"));
include(srcdir("helpers.jl"));
include(srcdir("stats.jl"));

# constructing kernel
N=2^12
r = Array(range(-pi,stop=pi,length=N),); 
eta = 0.0;
ξ = 2/3
expker = CovarianceKernel(r, eta, CovarianceCorrelation(expkernel), ξ, false); 
pwker = CovarianceKernel(r, eta, CovarianceCorrelation(piecewisekernel), ξ, false);
#expker = CovarianceKernel(r, eta, ξ);
logker = SingularKernel(r, eta);

# parseval and average
M = 10000
e1 = zeros(M)
e2 = zeros(M)
A = zeros(M)

for i in 1:M
    g1 = UnitaryWhiteNoise(div(N, 2) + 1);
    g2 = UnitaryWhiteNoise(div(N, 2) + 1);
    α = 0.6
    Γ = GmcNoise(logker,g1,g2,α);
    w = realization(pwker, Γ);
    
    e1[i] = 2*sum(abs.(Γ.wk[2:end] .* pwker.Lk[2:end]).^2) .+ abs.(Γ.wk[1] .* pwker.Lk[1]).^2
    e2[i] = mean(abs.(w).^2)
    
    A[i] = mean(w)
end
mean(A)
mean(e1)
mean(e2)

# 
α = 0.25
meanenergy(expker, 3000)
meanenergy(expker, 3000, α)



