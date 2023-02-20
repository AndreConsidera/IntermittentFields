using DrWatson
@quickactivate "IntermittentFields";
include(srcdir("IntermittentFields.jl"))
using .IntermittentFields
using FFTW
using Alert

# sizefunc
R12R13(p::AbstractArray{Float64}) = @. sqrt((p[1,:]-p[2,:])^2 + (p[1,:]-p[3,:])^2)

# all parameters
allparams = Dict(
    "np" => 3, 
    "ncor" => 400000,
    "nindep" => 1,     
    "α" => [0.0, 0.6],
    "ξ" => 2/3,
    "kernel" => piecewisekernel,
    "tmax" => 10,
    "N" => [2^i for i in 7:10],
    "beta" => 1.0,
    "planning_effort" => FFTW.MEASURE,
    "sizefunc" => R12R13,
    "dt" =>  1/(4 * 16 * 10^2),
    
)
dicts = dict_list(allparams);

function makesim(d::Dict)
    @unpack np, ncor, nindep, α, ξ, kernel, tmax, N, beta, planning_effort, sizefunc, dt = d
    t, p = dispersion(nindep, dt, N, α, ξ, tmax, beta, planning_effort, kernel, np, ncor, sizefunc)
    fulld = copy(d)
    fulld["t"] = t
    fulld["p"] = p
    return fulld
end

@alert for (i, d) in enumerate(dicts) 
    println("sim number=$i")
    fulldict = makesim(d)
    @tagsave(datadir("sims", "zeromodes", string(fulldict["kernel"]), "dispersion", savename("ET3l2normR12R13", d, "jld2")), fulldict, safe = DrWatson.readenv("DRWATSON_SAFESAVE", true))
end



#firstsim = readdir(datadir("sims","zeromodes","piecewisekernel", "dispersion"), join = true)[4]
#a=wload(firstsim)

#using DataFrames
#df = collect_results(datadir("sims","zeromodes","piecewisekernel", "dispersion"))


