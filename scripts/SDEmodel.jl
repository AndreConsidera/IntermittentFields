using DrWatson
@quickactivate "IntermittentFields";
include(srcdir("IntermittentFields.jl"))
using .IntermittentFields
using LinearAlgebra
using Statistics
using Distributions
using Interpolations
using Plots;plotly(size=(800,650))

function gffker(x, y, η)
    r=x-y
    @.(log(1/normevans(r, η))*(r <= 1))
end

function piecewisekernel(x, y, η, ξ)
    r=x-y
    @.((1. - normevans(r,η)^ξ) * (r<=1)) 
end

eta = 1e-12 
α = 0.3
ξ = 2/3
kappa=eta


tmax=500
np=2
nsamples = 1
et=ones(nsamples)*tmax
dt=1e-4
maxiter=Int(ceil(tmax/dt))
t= zeros(maxiter)
t= [i*dt for i in 0:maxiter-1] 
#p=zeros(maxiter, np)

Z=exp(α^2*log(2/eta))
for k in 1:nsamples  
    global p=zeros(maxiter, np)

    p[1,1] = 0
    p[1,2] = 2*eta
    exittime = 0 
    dw_til_1 = rand(Normal(0,1))*dt^0.5
    dw_til_2 = rand(Normal(0,1))*dt^0.5
    for i in 1:maxiter-1
        exittime = exittime + dt
        if abs(p[i,1] - p[i,2])>=1 
            break 
        end
        
        #singular kernel
        K_sq = [gffker(x, y, eta) for x in p[i,:], y in p[i,:]]
        if isposdef(K_sq)
            K, U = cholesky(K_sq)
        else
            break
        end
        #piecewise kernel
        L_sq = [piecewisekernel(x, y, eta, ξ) for x in p[i,:], y in p[i,:]]
        if isposdef(L_sq)
            L, U = cholesky(L_sq)
        else
            break
        end

        dw_1= rand(Normal(0,1))*dt^0.5
        dw_2=rand(Normal(0,1))*dt^0.5
        
        k_noise_1 = rand(Normal(0,1))*(2*kappa*dt)^0.5
        k_noise_2 = rand(Normal(0,1))*(2*kappa*dt)^0.5
        p[i+1, 1] = p[i, 1] + L[1,1]*(1/Z)*exp(α*(K[1,1]*dw_til_1+K[1,2]*dw_til_2))*dw_1 + L[1,2]*(1/Z)*exp(α*(K[2,1]*dw_til_1+K[2,2]*dw_til_2))*dw_2 + k_noise_1
        p[i+1, 2] = p[i, 2] + L[2,1]*(1/Z)*exp(α*(K[1,1]*dw_til_1+K[1,2]*dw_til_2))*dw_1 + L[2,2]*(1/Z)*exp(α*(K[2,1]*dw_til_1+K[2,2]*dw_til_2))*dw_2 + k_noise_2

    end
    et[k] = exittime
end

mean(et)
p1 = plot(t[t.<et[end]-dt],[p[t.<et[end]-dt,1],p[t.<et[end]-dt,2]],label="");
p2 = plot(t[t.<et[end]-dt],abs.(p[t.<et[end]-dt,1] .- p[t.<et[end]-dt,2]),label="");
layout = @layout [a ; c]
plot(p1, p2, layout = layout)

α_c= sqrt((1-ξ)/2)
D0 = 1
αmean(λ, ξ, α) = (λ^2)/(2 * ((2/(2-(ξ+2*α^2)))*(1-(ξ+2*α^2)^2)/(1+(ξ+2*α^2))) * (D0*(1-((ξ+2*α^2)/2))^2))

#δ = circshift(p[:,1],-1) - p[:,1]
#plot(t[1:end-1],δ[1:end-1])
#xinterp = minimum(x):0.001:maximum(x)
#p3 = plot(xinterp, gmcinterp(xinterp));

