using Distributed
procs = addprocs(8)
@everywhere using DrWatson
@everywhere @quickactivate "IntermittentFields";
@everywhere include(srcdir("IntermittentFields.jl"))
@everywhere using .IntermittentFields
@everywhere using FFTW
using Alert

Np = Int(ceil(8*13/length(procs)))
dt = 1/(4 * 16 * 10^2)
N = 2^8
α = 0.2
ξ = 2/3
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

nproc = length(procs)
params = @strdict Np dt N α ξ tmax beta nproc;
out = @strdict t t2 t4 t8 t16 t32 t64 p p2 p4 p8 p16 p32 p64;

a2 = @spawnat procs[1] doexittime2(Np, dt, N, α, ξ, tmax, beta, planning_effort)
a3 = @spawnat procs[2] doexittime2(Np, dt, N, α, ξ, tmax, beta, planning_effort)
a4 = @spawnat procs[3] doexittime2(Np, dt, N, α, ξ, tmax, beta, planning_effort)
a5 = @spawnat procs[4] doexittime2(Np, dt, N, α, ξ, tmax, beta, planning_effort)
a6 = @spawnat procs[5] doexittime2(Np, dt, N, α, ξ, tmax, beta, planning_effort)
a7 = @spawnat procs[6] doexittime2(Np, dt, N, α, ξ, tmax, beta, planning_effort)
a8 = @spawnat procs[7] doexittime2(Np, dt, N, α, ξ, tmax, beta, planning_effort)
a9 = @spawnat procs[8] doexittime2(Np, dt, N, α, ξ, tmax, beta, planning_effort)

path=datadir("sims", "2particles","non_rescaled_kernel", savename("ET2",params,"jld2"));
@alert begin
    t[1,:],t2[1,:],t4[1,:],t8[1,:],t16[1,:],t32[1,:],t64[1,:], p[1, :, :], p2[1, :, :],p4[1, :, :],p8[1, :, :],p16[1, :, :],p32[1, :, :],p64[1, :, :] = fetch(a2)
    t[2,:],t2[2,:],t4[2,:],t8[2,:],t16[2,:],t32[2,:],t64[2,:], p[2, :, :], p2[2, :, :],p4[2, :, :],p8[2, :, :],p16[2, :, :],p32[2, :, :],p64[2, :, :] = fetch(a3)
    t[3,:],t2[3,:],t4[3,:],t8[3,:],t16[3,:],t32[3,:],t64[3,:], p[3, :, :], p2[3, :, :],p4[3, :, :],p8[3, :, :],p16[3, :, :],p32[3, :, :],p64[3, :, :] = fetch(a4)
    t[4,:],t2[4,:],t4[4,:],t8[4,:],t16[4,:],t32[4,:],t64[4,:], p[4, :, :], p2[4, :, :],p4[4, :, :],p8[4, :, :],p16[4, :, :],p32[4, :, :],p64[4, :, :] = fetch(a5)
    t[5,:],t2[5,:],t4[5,:],t8[5,:],t16[5,:],t32[5,:],t64[5,:], p[5, :, :], p2[5, :, :],p4[5, :, :],p8[5, :, :],p16[5, :, :],p32[5, :, :],p64[5, :, :] = fetch(a6)
    t[6,:],t2[6,:],t4[6,:],t8[6,:],t16[6,:],t32[6,:],t64[6,:], p[6, :, :], p2[6, :, :],p4[6, :, :],p8[6, :, :],p16[6, :, :],p32[6, :, :],p64[6, :, :] = fetch(a7)
    t[7,:],t2[7,:],t4[7,:],t8[7,:],t16[7,:],t32[7,:],t64[7,:], p[7, :, :], p2[7, :, :],p4[7, :, :],p8[7, :, :],p16[7, :, :],p32[7, :, :],p64[7, :, :] = fetch(a8)
    t[8,:],t2[8,:],t4[8,:],t8[8,:],t16[8,:],t32[8,:],t64[8,:], p[8, :, :], p2[8, :, :],p4[8, :, :],p8[8, :, :],p16[8, :, :],p32[8, :, :],p64[8, :, :] = fetch(a9)
    
    safesave(path, out)    
end


#=
using JLD2
f = jldopen(path, "r")
p=f["p4"]
@show p
close(f)
mean(A)
=#