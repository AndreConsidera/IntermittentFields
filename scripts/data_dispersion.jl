using DrWatson
@quickactivate "IntermittentFields";
include(srcdir("IntermittentFields.jl"))
using .IntermittentFields
using FFTW
using Alert

using Statistics
using StatsPlots;plotly()
using StatsBase
using LinearAlgebra

np = 3     # 2 -> pairs , 3 -> triplets
sizefunc(p::AbstractArray{Float64}) = @. sqrt((p[1,:]-p[2,:])^2 + (p[1,:]-p[3,:])^2)
#sizefunc(p::AbstractArray{Float64}) = @. sqrt((p[1,:]-p[2,:])^2)
ncor = 400000    # Number of correlated samples
nindep = 1   # Number of independent samples
dt = 1/(4 * 16 * 10^2)
N = 2^8
α = 0.0
ξ = 2/3
kernel = piecewisekernel
tmax = 10.
beta = 1.
planning_effort = FFTW.MEASURE

params = @strdict np dt N α ξ tmax beta ncor nindep;
out = @strdict t p ;
path = datadir("sims", "zeromodes", string(kernel), "dispersion", savename("ET3l2normR12R13", params, "jld2"));
 
t, p = dispersion(nindep, dt, N, α, ξ, tmax, beta, planning_effort, kernel, np, ncor, sizefunc)
safesave(path, out)    
