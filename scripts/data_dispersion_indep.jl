using DrWatson
@quickactivate "IntermittentFields";
include(srcdir("IntermittentFields.jl"))
using .IntermittentFields
using FFTW
using Alert
using Dates


et = 10*ones(100)
et2 = copy(et)
# condition return vector of Bool
function condition1(p, t)
    [sqrt((p[1,i]-p[2,i])^2 + (p[1, i]-p[3, i])^2) >= 1 for i in 1:100]
end

function condition2(p, t)
    [sqrt((p[1,i]-p[2,i])^2 + (p[1, i]-p[3, i])^2) >= 1/2 for i in 1:100]
end


function affect1(p, t, idx)
    et[idx] = min(et[idx], t[idx])
end

function affect2(p, t, idx)
    et2[idx] = min(et2[idx], t[idx])
end

cb1 = VectorCallback(condition1, affect1)
cb2 = VectorCallback(condition2, affect2)

cs = VectorCallbackSet((cb1,cb2))


#======================== callbacks for ncor = 1
exittime = 10*ones(10)
exittime2= 10*ones(10)

# callback1 ==============================
condition1(p, t, k) = @. sqrt((p[1,1]-p[2,1])^2 + (p[1,1]-p[3,1])^2) >= 1
function affect1(p, t, k)
    exittime[k] = min(exittime[k], t[1])
end
cb1 = DiscreteCallback(condition1, affect1)
# ===================================

# callback2 ==============================
condition2(p, t, k) = @. sqrt((p[1,1]-p[2,1])^2 + (p[1,1]-p[3,1])^2) >= 1/64
function affect2(p, t, k)
    exittime2[k] = min(exittime2[k], t[1])
end
cb2 = DiscreteCallback(condition2, affect2)
# ===================================
    
cs = CallbackSet((cb1,cb2))
=============================================================#

# sizefunc
R12R13(p::AbstractArray{Float64}) = @. sqrt((p[1,:]-p[2,:])^2 + (p[1,:]-p[3,:])^2)
    
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

# all parameters
allparams = Dict(
    "np" => 3, 
    "ncor" => 100,
    "nindep" => 1,     
    "α" => [0.6],
    "ξ" => 2/3,
    "kernel" => piecewisekernel,
    "tmax" => 10,
    "N" => [2^i for i in 7:7],
    "beta" => 1.0,
    "planning_effort" => FFTW.MEASURE,
    "sizefunc" => R12R13,
    "dt" =>  1/(4 * 16 * 10^2),
    
)
dicts = dict_list(allparams);

function makesim(d::Dict)
    @unpack np, ncor, nindep, α, ξ, kernel, tmax, N, beta, planning_effort, sizefunc, dt = d
    t, p = dispersion(nindep, dt, N, α, ξ, tmax, beta, planning_effort, kernel, np, ncor, sizefunc, logmessage, cs)
    fulld = copy(d)
    fulld["t"] = t
    fulld["p"] = p
    return fulld
end

@alert for (i, d) in enumerate(dicts) 
    println("sim number=$i/$(length(dicts))")
    fulldict = makesim(d)
    #@tagsave(datadir("sims", "zeromodes", string(fulldict["kernel"]), "dispersion", savename("ET3l2normR12R13", d, "jld2")), fulldict, safe = DrWatson.readenv("DRWATSON_SAFESAVE", true))
end


sim = filter(x->occursin("ncor=40_", x), readdir(datadir("sims","zeromodes","piecewisekernel","dispersion"), join = true))
a=wload(sim[1])
t=vec(a["t"])
et
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