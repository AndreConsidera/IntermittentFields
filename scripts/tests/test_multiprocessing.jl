using DrWatson;
@quickactivate("IntermittentFields");
include(srcdir("intermittentfields_mod.jl"))
import FFTW
using .IntermittentFields

Np = 320
dt = 1/200
N = 2^12
α = 0.2
ξ = 2/3
tmax = 10 
beta = 1.
planning_effort = FFTW.MEASURE
@time doexittime(Np, dt, N, α, ξ, tmax, beta)

# mutiprocess 8 cores

using Distributed
procs = addprocs(8)
@everywhere using DrWatson
@everywhere @quickactivate("IntermittentFields");
@everywhere include(srcdir("intermittentfields_mod.jl"))
@everywhere using .IntermittentFields
@everywhere using FFTW

ET = zeros(8, Int(Np/8));
@time begin
    a2 = @spawnat procs[1] doexittime(Int(ceil(Np/8)), dt, N, α, ξ, tmax, beta)
    a3 = @spawnat procs[2] doexittime(Int(ceil(Np/8)), dt, N, α, ξ, tmax, beta)
    a4 = @spawnat procs[3] doexittime(Int(ceil(Np/8)), dt, N, α, ξ, tmax, beta)
    a5 = @spawnat procs[4] doexittime(Int(ceil(Np/8)), dt, N, α, ξ, tmax, beta)
    a6 = @spawnat procs[5] doexittime(Int(ceil(Np/8)), dt, N, α, ξ, tmax, beta)
    a7 = @spawnat procs[6] doexittime(Int(ceil(Np/8)), dt, N, α, ξ, tmax, beta)
    a8 = @spawnat procs[7] doexittime(Int(ceil(Np/8)), dt, N, α, ξ, tmax, beta)
    a9 = @spawnat procs[8] doexittime(Int(ceil(Np/8)), dt, N, α, ξ, tmax, beta)

    ET[1,:]=fetch(a2)
    ET[2,:]=fetch(a3)
    ET[3,:]=fetch(a4)
    ET[4,:]=fetch(a5)
    ET[5,:]=fetch(a6)
    ET[6,:]=fetch(a7)
    ET[7,:]=fetch(a8)
    ET[8,:]=fetch(a9)  
end