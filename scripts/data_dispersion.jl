using Distributed
nproc = 8
procs = addprocs(nproc)
@everywhere using DrWatson
@everywhere @quickactivate "IntermittentFields";
@everywhere include(srcdir("IntermittentFields.jl"))
@everywhere using .IntermittentFields
@everywhere using FFTW


# sizefunc
@everywhere R12R13(p::AbstractArray{<:Real}) = @. sqrt((p[1,:]-p[2,:])^2 + (p[1,:]-p[3,:])^2)
@everywhere R12(p::AbstractArray{<:Real}) = @. sqrt((p[1,:]-p[2,:])^2)


# all parameters
d = Dict(
    "np" => [2], 
    "ncor" => [100],
    "nindep" => [1],
    "α" => [0.3],
    "ξ" => [2/3],
    "kernel" => piecewisekernel,
    "tmax" => 10,
    "λ" => [1/2^0],
    "N" => [2^i for i in 12:12],
    "dt"=> 1e-3,
    "beta" => 1.0,
    "planning_effort" => FFTW.MEASURE,
    "sizefunc" => [@onlyif("np" == 2, R12), @onlyif("np" == 3, R12R13)], 
)

#d["dt"] = vec([((4 * pi ./ d["N"][i])./d["λ"][j]).^(2 .-d["ξ"][k]) for i in eachindex(d["N"]), j in eachindex(d["λ"]), k in eachindex(d["ξ"]) ])

dicts = dict_list(d);

function makesim(d::Dict)
    @unpack np, ncor, nindep, α, ξ, kernel, tmax, λ, N, dt, beta, planning_effort, sizefunc = d
    t = zeros(ncor, nindep, nproc);
    p = zeros(np, ncor, nindep, nproc);
    a = Array{Future}(undef, nproc)
    for i in 1:nproc
        a[i] = @spawnat procs[i]  dispersion(nindep, dt, N, α, ξ, tmax, beta, planning_effort, kernel, np, ncor, sizefunc, λ)
    end
    for i in 1:nproc
        t[:,:,i], p[:, :, :, i] = fetch(a[i])
    end
    fulld = copy(d)
    fulld["t"] = t
    fulld["p"] = p
    fulld["nproc"] = nproc
    return fulld
end

for (i, d) in enumerate(dicts) 
    println("sim number=$i/$(length(dicts))")
    println("parameters:", savename(d),"_nproc=$nproc")
    global fulldict = makesim(d)
#    @tagsave(datadir("sims", "dispersion", "2particle_scaling", savename("noise=scaling", d, "jld2")), fulldict, safe = DrWatson.readenv("DRWATSON_SAFESAVE", true))
end
