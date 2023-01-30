using DrWatson
@quickactivate("IntermittentFields");
include(srcdir("IntermittentFields.jl"))
using .IntermittentFields
using Random
using Statistics
using Distributions
using Plots;plotly()
using Interpolations
using JLD2
using LaTeXStrings
using FFTW
using Polynomials

#plot kernels
N=2^14
r = Array(range(0,stop=L,length=N),); 
eta = 0.0 #4*pi/N;
ξ = 1/3

expker = CovarianceKernel(r, eta, ξ);
resexp = CovarianceKernel(r, eta, CovarianceCorrelation(expkernel), ξ, false); 
pwker = CovarianceKernel(r, eta, CovarianceCorrelation(piecewisekernel), ξ, false);
logker = SingularKernel(r, eta);

# plot kernels 
begin
    plot(ylims=(0,4),xlims=(0,1),title="Kernels",xlabel="r")
    plot!(r,expkernel(r,ξ), label="Exp")
    plot!(r,piecewisekernel(r,ξ), label="Piecewise")
#    plot!(logker.r,logker.cr, label="Log")
end
#savefig(p,plotsdir(savename("kernels",d,"png")));

c0=exp(eta/2)^(ξ)

begin
    plot(ylims = (0,1), xlims = (1e-4,1), title = "Kernels", xlabel = "r",xaxis=:log,yaxis=:log)

   # plot!(expker.r[r.>0], expker.r[r.>0].^(-ξ) .*(1 .-expker.cr[r.>0]), label="Exp",linealpha = 0.2, lw = 6,xaxis=:linear,yaxis=:linear)
   # plot!(expker.r, expker.Lr_sq, label=L"L*L")
    plot!(resexp.r[resexp.r.>0], resexp.r[resexp.r.>0].^(0).*(maximum(resexp.Lr_sq) .-resexp.Lr_sq[resexp.r.>0]), label=L"not resc")
   # plot!(pwker.r, pwker.cr)
    plot!(resexp.r[resexp.r.>0], resexp.r[resexp.r.>0] .^(0).*(maximum(resexp.cr).-resexp.cr[resexp.r.>0]),lw=15,alpha=0.4)
    plot!(resexp.r[resexp.r.>0], resexp.r[resexp.r.>0] .^(ξ),lw=3,alpha=1)

    # plot!(pwker.r[r.>0],exp.(-1.36*pwker.r[r.>0] .^ξ))
    #plot!(pwker.r, expker.cr, label="piecewise",linealpha = 0.2, lw = 6)
    #plot!(pwker.r, expker.Lr_sq, label=L"PW L*L")
end

# Gaussian fields
g1 = UnitaryWhiteNoise(div(N, 2) + 1);
u = realization(expker, g1);
plot(r, u)

#Intermittent fields
g2 = UnitaryWhiteNoise(div(N, 2) + 1);
α = 0.25
Γ = GmcNoise(logker,g1,g2,α);
v = realization(expker, Γ);
plot(r,[u, v])

# plot sp
p = 4
m = 1000
α = 0.2
sp_abs, sp, l = structurefunc(expker, m, p, α);
sp_abs, sp, l = structurefunc(expker_rescaled, m, p, α);
plot(l, sp_abs,xaxis=:log, yaxis=:log, markershape=:circle );
tmp = linear_interpolation(l,sp_abs)(1);
plot!(l, tmp .* l.^(p*expker_rescaled.ξ/2))

# PDF of δv
δ = δv(expker_rescaled, 1000,0.25);
P = histogram(yaxis=:log,ylims=(1e-1,1e6),xlims=(-7,7));
for i in axes(δ,1)[1:2:end]
    histogram!((δ[i,:]))
end
display(P)

# ζ(p) ∼ lognormal
function ζ(p, m, α)
    sp_abs, sp, l = structurefunc(expker_rescaled, m, p, α);

    pol=Polynomials.fit(log.(l), log.(sp_abs), 1)
    coeffs(pol)[2]
    return coeffs(pol)[2]
end

α=0.25
z=[ζ(p, 1000, α) for p in 1:8]
x=[i for i in 1:8]
scatter(1:8,z,markershape=:circle);
plot!(x,x.*ξ/2);
plot!(x,x.*ξ/2 - x.*(x./2 .-1).*α^2)

