using DrWatson
using Random
using Statistics
using Distributions
using Plots
using DifferentialEquations
using Interpolations
using JLD2
using Test
using LaTeXStrings
using Distributed
using ProgressMeter

@quickactivate("Intermittency Paradox");
include(srcdir("general_lib.jl"));
include(srcdir("compute_kernels.jl"));

# change distribution of particle positions

N=2^12
d = Dict(
    :r =>Array(range(-pi,stop=pi,length=N)), 
    :eta=>4*pi/N,
    :ker=>ExpKernel(ξ=0.6),
    :rescale=>false
    );

#constructing field    
cr,ck,Lr,Lr_sq,Lk=compute_kernels(;d...);L=CovarianceConvolutionKernel(Lr,Lr_sq,Lk);f=GaussianField(L);

# drift and sigma
drift(u,p,t) = 0
function sig(u,p,t) 
    wk=(1/2^0.5)*(rand(Normal(0,1),length(d[:r])÷2 + 1) +im*rand(Normal(0,1),length(d[:r])÷2 + 1));
    noise=WhiteNoise(wk)
    v=realization(f,noise)
    uinterp = linear_interpolation(d[:r],v,extrapolation_bc=Periodic())
    
    if p[:savefields]==true
        jldopen(datadir("savedfields.jld2"), "a+") do f
            if isempty(f)==true
                f["field_1"]=uinterp
            else
                fieldcounter=split(keys(f)[end],"_")[2]
                f["field_"*string(parse(Int,fieldcounter)+1)] = uinterp
            end
        end  
    end
    # kappa for point splitting

    return uinterp.(u) .+ sqrt(2*p[:kappa])*rand(Normal(0,1),size(u,1))
end

# parameters
p=Dict(:kappa=>d[:eta]^1,:savefields=>false)
Np=Int(16e2);u0=[0.,2. *d[:eta]];
dt = 1e-2/16.;tmax=10.;tspan = (0.0,tmax);
# Noise
Wtmp=[(i*1.)*dt^0.5 for i in eachindex(0:dt:tmax)];W=NoiseGrid(0:dt:tmax,Wtmp);
# solve SDEProblem
prob = SDEProblem(drift,sig,u0,tspan,p,noise=W);

times=zeros(Np)
M=zeros(5)
for k in eachindex(M)
    @showprogress 1 "Computing...$k"    for i in eachindex(times)
        sol = solve(prob,EM(false),dt=dt);
        asol=convert(Array,sol);
        rr=abs.(asol[1,:]-asol[2,:]);
        if length(sol.t[rr.>1])==0
            times[i]=tmax
        else
            times[i]=sol.t[rr.>1][1]   
        end   
    end
    M[k]=mean(times)
end

std(M)


# test method
@testset "euler method via SDE Solver" for eq in (isequal,) 
    file=jldopen(datadir("savedfields.jld2"),"r")
    for i in (1:size(sol.u)[1]-1)
        v=read(file,keys(file)[i])
        @test eq(sol.u[i+1] , sol.u[i] + v(sol.u[i])*dt^0.5)
    end
    close(file)
end

