using DrWatson
@quickactivate "IntermittentFields";
include(srcdir("IntermittentFields.jl"))
using .IntermittentFields
using LinearAlgebra
using Statistics
using Plots;plotly(size=(800,650))
using DataFrames
using Statistics
using StatsBase
using Distributions
using Interpolations
import Random

N=2^14
tmax=10
#dt=1e-2
nindep = 5000
α = 0.2

dts = [1e-2]
m=zeros(length(dts))
for (j,dt) in enumerate(dts)

#Random.seed!(1234567)
 

dx = 2 * pi/N    
eta = 2 * dx
r0 = 2 * eta
r = [-pi + i * dx for i in 0:N-1]
kappa = eta
#instantiating kernels and FFT
logker = SingularKernel(r, eta)
push!(r, pi)
maxiter = Int(ceil(tmax/dt))
exittime = zeros(nindep);
global t =[i*dt for i in 0:maxiter-1]

finalp = zeros(nindep);
global traj = zeros(maxiter,nindep);
p = zeros(maxiter)
#p[2] = p[1] + r0



for k in 1:nindep
    #setting environment
    g2 = UnitaryWhiteNoise(div(length(logker.r), 2) + 1)
    σ_sq = 2 * sum(abs.(logker.Lk[2:end]).^2)
    gmc = exp.(α .* realization(logker, g2) .- α^2 .* σ_sq)
    gmcinterp = linear_interpolation(r, push!(gmc, gmc[1]), extrapolation_bc = Periodic())

    p[1] = 3 #2*rand(Uniform(0, 1))
    for i in 1:maxiter-1
        #t[i+i] = t[i] + dt
        dw_1= rand(Normal(0,1))*dt^0.5
        p[i+1] = p[i] + gmcinterp(p[i])*dw_1 

    end
    #exittime[k] = t
    finalp[k] = p[end]
    traj[:,k]=p
end



m[j]=mean(abs.(finalp))
end

plot(t,[mean(traj,dims=2)])
#plot!(t,traj[:,4])

scatter(dts,m,xaxis=:log,yaxis=:log,ylims=(1e-3,1))

plot(r,gmcinterp.(r))
plot!(t,p)

# OU ---