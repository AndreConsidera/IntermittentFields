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
input = Dict(
    "np" => [3], 
    "ncor" => [500],
    "nindep" => [100],
    "α" => [0.0],
    "ξ" => [2/3],
    "kernel" => piecewisekernel,
    "tmax" => 10,
    "λ" => [1],
    "N" => [2^i for i in 7:9],
    "beta" => 1.0,
    "planning_effort" => FFTW.MEASURE,
    "sizefunc" => [@onlyif("np" == 2, R12), @onlyif("np" == 3, R12R13)],
    "dt" =>  1/(4* 16 * 10^2), 
)
dicts = dict_list(input);

function makesim(d::Dict)
    @unpack np, ncor, nindep, α, ξ, kernel, tmax, λ, N, beta, planning_effort, sizefunc, dt = d
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
    @tagsave(datadir("sims", "zeromodes", string(fulldict["kernel"]), "dispersion", savename("ET3l2normR12R13", d, "jld2")), fulldict, safe = DrWatson.readenv("DRWATSON_SAFESAVE", true))
end
