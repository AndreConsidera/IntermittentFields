using DrWatson
@quickactivate "IntermittentFields";
include(srcdir("IntermittentFields.jl"))
using .IntermittentFields
using FFTW
using Alert
using Dates

# sizefunc
R12R13(p::AbstractArray{Real}) = @. sqrt((p[1,:]-p[2,:])^2 + (p[1,:]-p[3,:])^2)
R12(p::AbstractArray{Real}) = @. sqrt((p[1,:]-p[2,:])^2)

# all parameters
input = Dict(
    "np" => [2, 3], 
    "ncor" => 100,
    "nindep" => 10,     
    "α" => [0.0, 0.6],
    "ξ" => [1/3, 2/3],
    "kernel" => piecewisekernel,
    "tmax" => 10,
    "λ" => 1,
    "N" => [2^i for i in 7:9],
    "beta" => 1.0,
    "planning_effort" => FFTW.MEASURE,
    "sizefunc" => [@onlyif("np" == 2, R12), @onlyif("np" == 3, R12R13)],
    "dt" =>  1/(4* 16 * 10^2), 
)
dicts = dict_list(input);

# logmessage
function logmessage(i, maxiter, particles_stuck) 
    # current time
    time = Dates.format(now(UTC), dateformat"yyyy-mm-dd HH:MM:SS")
    # memory the process is using 
    maxrss = "$(round(Sys.maxrss()/1048576, digits=2)) MiB"
    percent = (i/maxiter)*100
    logdata = (; 
    percent, # iteration i
    particles_stuck, 
    maxrss) # lastly the amount of memory being used
    
    if mod(i, maxiter/10)==0
        println(savename(time, logdata; connector = " | ", equals = " = ", sort = false, digits = 3))
    end
    
end


function makesim(d::Dict)
    @unpack np, ncor, nindep, α, ξ, kernel, tmax, λ, N, beta, planning_effort, sizefunc, dt = d
    t, p = dispersion(nindep, dt, N, α, ξ, tmax, beta, planning_effort, kernel, np, ncor, sizefunc, logmessage, λ)
    fulld = copy(d)
    fulld["t"] = t
    fulld["p"] = p
    return fulld
end

@alert for (i, d) in enumerate(dicts) 
    println("sim number=$i/$(length(dicts))")
    global fulldict = makesim(d)
    @tagsave(datadir("sims", "zeromodes", string(fulldict["kernel"]), "dispersion", savename("ET3l2normR12R13", d, "jld2")), fulldict, safe = DrWatson.readenv("DRWATSON_SAFESAVE", true))
end



#fulldict["t"]
#sum(fulldict["et"][2,:] .< vec(fulldict["t"])) == input["ncor"]


#sim = filter(x->occursin("ncor=40_", x), readdir(datadir("sims","zeromodes","piecewisekernel","dispersion"), join = true))
#a=wload(sim[1])
#t=vec(a["t"])
#et
#firstsim = readdir(datadir("sims","zeromodes","piecewisekernel", "dispersion"), join = true)

#=
using DataFrames
using Statistics

meanexittime = [:mean_exittime => d -> mean(d["t"]), :shape => d -> mean(abs.(d["p"][1,:,:] - d["p"][2,:,:]))]

bl = ["gitcommit","tmax","p", "planning_effort", "t","sizefunc","script"]
df = collect_results(datadir("sims","zeromodes","piecewisekernel", "dispersion"); black_list = bl, special_list = meanexittime)

vscodedisplay(df, "results_piecewise")
scatter(df.N, df.shape)
=#