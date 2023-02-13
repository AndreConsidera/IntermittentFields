using Distributed
nproc = 8
procs = addprocs(nproc)
@everywhere using DrWatson
@everywhere @quickactivate "IntermittentFields";
@everywhere include(srcdir("IntermittentFields.jl"))
@everywhere using .IntermittentFields
@everywhere using FFTW
using Alert

Np = Int(ceil(nproc * 9/length(procs)))
dt = 1/(4 * 16 * 10^2)
N = 2^7
α = 0.0
ξ = 2/3
kernel = piecewisekernel
tmax = 10.
beta = 1.
planning_effort = FFTW.MEASURE

t = zeros(length(procs), Np);
t2 = zeros(length(procs), Np);
t4 = zeros(length(procs), Np);
t8 = zeros(length(procs), Np);
t16 = zeros(length(procs), Np);
t32 = zeros(length(procs), Np);
t64 = zeros(length(procs), Np);

p = zeros(length(procs), 2, Np);
p2 = zeros(length(procs), 2, Np);
p4 = zeros(length(procs), 2, Np);
p8 = zeros(length(procs), 2, Np);
p16 = zeros(length(procs), 2, Np);
p32 = zeros(length(procs), 2, Np);
p64 = zeros(length(procs), 2, Np);

# input output dict
params = @strdict Np dt N α ξ tmax beta nproc;
out = @strdict t t2 t4 t8 t16 t32 t64 p p2 p4 p8 p16 p32 p64;

a = Array{Future}(undef, nproc)
for i in 1:nproc
    a[i] = @spawnat procs[i]  doexittime2(Np, dt, N, α, ξ, tmax, beta, planning_effort, kernel)
end

path = datadir("sims", "2particles", string(kernel), savename("ET2", params, "jld2"));
@alert begin    
    for i in 1:nproc
        t[i,:],t2[i,:],t4[i,:],t8[i,:],t16[i,:],t32[i,:],t64[i,:], p[i, :, :], p2[i, :, :],p4[i, :, :],p8[i, :, :],p16[i, :, :],p32[i, :, :],p64[i, :, :] = fetch(a[i])
    end
    safesave(path, out)
end

#=
using JLD2
f = jldopen(path, "r")
A=f["t2"]
@show A
close(f)
mean(A)
=#
