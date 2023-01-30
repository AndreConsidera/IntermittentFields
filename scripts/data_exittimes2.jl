using Distributed
procs = addprocs(8)
@everywhere using DrWatson
@everywhere @quickactivate "Intermittency Paradox";
@everywhere include(srcdir("IntermittentFields.jl"))
@everywhere using .IntermittentFields
@everywhere using FFTW
using Alert

Np = Int(ceil(8*200/length(procs)))
dt = 1/(4 * 16 * 10^2)
N = 2^12
α = 0.0
ξ = 1/3
tmax = 10.
beta = 1.
planning_effort = FFTW.MEASURE

ET = zeros(length(procs),Np);
nproc = length(procs)
params = @strdict Np dt N α ξ tmax beta nproc;
out = @strdict ET;

a2 = @spawnat procs[1] doexittime(Np, dt, N, α, ξ, tmax, beta, planning_effort)
a3 = @spawnat procs[2] doexittime(Np, dt, N, α, ξ, tmax, beta, planning_effort)
a4 = @spawnat procs[3] doexittime(Np, dt, N, α, ξ, tmax, beta, planning_effort)
a5 = @spawnat procs[4] doexittime(Np, dt, N, α, ξ, tmax, beta, planning_effort)
a6 = @spawnat procs[5] doexittime(Np, dt, N, α, ξ, tmax, beta, planning_effort)
a7 = @spawnat procs[6] doexittime(Np, dt, N, α, ξ, tmax, beta, planning_effort)
a8 = @spawnat procs[7] doexittime(Np, dt, N, α, ξ, tmax, beta, planning_effort)
a9 = @spawnat procs[8] doexittime(Np, dt, N, α, ξ, tmax, beta, planning_effort)

path=datadir("sims", "2particles", savename("ET2",params,"jld2"));
@alert begin
    ET[1,:]=fetch(a2)
    ET[2,:]=fetch(a3)
    ET[3,:]=fetch(a4)
    ET[4,:]=fetch(a5)
    ET[5,:]=fetch(a6)
    ET[6,:]=fetch(a7)
    ET[7,:]=fetch(a8)
    ET[8,:]=fetch(a9)    
    
    safesave(path, out)    
end

#=
f = jldopen(path, "r")
@show ET 
close(f)
A=f["ET"]
mean(A)
=#