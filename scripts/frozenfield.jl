using Distributed
nproc = 8
procs = addprocs(nproc)
@everywhere using DrWatson
@everywhere @quickactivate "IntermittentFields";
@everywhere include(srcdir("IntermittentFields.jl"))
@everywhere using .IntermittentFields
@everywhere using FFTW
@everywhere using Distributions
@everywhere using Interpolations
@everywhere using LinearAlgebra

@everywhere function piecewisekernel(x, y, η, ξ)
    r=x-y
    @.((1. - normevans(r,η)^ξ) * (r<=1)) 
end

@everywhere function expkernel(x, y, η, ξ)
    r=x-y
    @.(exp(-normevans(r,η)^ξ))
end


@everywhere function hybridgmc(nindep::Integer,dt::Real, N::Real ,α::Real, ξ::Real, tmax::Real, beta::Real, np::Integer)
    dx = 2 * pi/N    
    eta = 2 * dx
    r0 = 2 * eta
    r = [-pi + i * dx for i in 0:N-1]
    kappa = eta^beta
    #instantiating kernels and FFT
    logker = SingularKernel(r, eta)
    push!(r, pi)
    maxiter = ceil(tmax/dt);
    exittime = zeros(nindep);
    finalp = zeros(np, nindep);
    for k in 1:nindep 
        p = zeros(np)
        p[1] = rand(Uniform(-pi, pi))
        p[2] = p[1] + r0
        t = 0 
        
        #setting environment
        g2 = UnitaryWhiteNoise(div(length(logker.r), 2) + 1)
        σ_sq = 2 * sum(abs.(logker.Lk[2:end]).^2)
        gmc = exp.(α .* realization(logker, g2) .- α^2 .* σ_sq)
        gmcinterp = linear_interpolation(r, push!(gmc, gmc[1]), extrapolation_bc = Periodic())

        for i in 1:maxiter-1
            if abs(p[1] - p[2])>=1 
                break 
            end
            t = t + dt

            #piecewise kernel
            L_sq = [piecewisekernel(x, y, eta, ξ) for x in p, y in p]
            if isposdef(L_sq)
                L, U = cholesky(L_sq)
            else
                L = sqrt(L_sq)
            end

            dw_1= rand(Normal(0,1))*dt^0.5
            dw_2=rand(Normal(0,1))*dt^0.5
            k_noise_1 = rand(Normal(0,1))*(2*kappa*dt)^0.5
            k_noise_2 = rand(Normal(0,1))*(2*kappa*dt)^0.5

            p[1] = p[1] + L[1,1]*gmcinterp(p[1])*dw_1 + L[1,2]*gmcinterp(p[2])*dw_2 + k_noise_1
            p[2] = p[2] + L[2,1]*gmcinterp(p[1])*dw_1 + L[2,2]*gmcinterp(p[2])*dw_2 + k_noise_2

        end
        exittime[k] = t
        finalp[:, k] = p
    end
    return exittime, finalp
end

# all parameters
d = Dict{String, Any}(
"np" => 2, 
"nindep" => 10,
"α" => 0.2,
"ξ" => 2/3,
"tmax" => 128,
"N" => 2^14,
"dt"=> 1e-4,
"beta" => 1.0,
"placeholder"=>"togetcorrecttype "
)

dicts = dict_list(d);
    
function makesim(d::Dict)
    @unpack np, nindep, α, ξ, tmax, N, dt, beta = d
    t = zeros(nindep, nproc);
    p = zeros(np, nindep, nproc);
    a = Array{Future}(undef, nproc)
    for i in 1:nproc
        a[i] = @spawnat procs[i]  hybridgmc(nindep::Integer,dt::Real, N::Real ,α::Real, ξ::Real, tmax::Real, beta::Real, np::Integer)
    end
    for i in 1:nproc
        t[:,i], p[:, :, i] = fetch(a[i])
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
    #@tagsave(datadir("sims", "dispersion", "frozenfield", savename("frozen", d, "jld2")), fulldict, safe = DrWatson.readenv("DRWATSON_SAFESAVE", true))
end     



#
p1 = plot(t[t.<et[end]-dt],[p[t.<et[end]-dt,1],p[t.<et[end]-dt,2]],label="");
p2 = plot(t[t.<et[end]-dt],abs.(p[t.<et[end]-dt,1] .- p[t.<et[end]-dt,2]),label="");
layout = @layout [a ; c]
plot(p1, p2, layout = layout)
α_c= sqrt((1-ξ)/4)
D0 = 1
αmean(λ, ξ, α) = (λ^2)/(2 * ((2/(2-(ξ+4*α^2)))*(1-(ξ+4*α^2)^2)/(1+(ξ+4*α^2))) * (D0*(1-((ξ+4*α^2)/2))^2))
αmean(1,1/3,0.1)


# data analysis
using DrWatson
@quickactivate "IntermittentFields";
using LinearAlgebra
using Statistics
using Plots;plotly(size=(800,650))
using DataFrames
using Statistics
using StatsBase
white_list = ["N", "α", "p", "t", "ξ", "nindep","np", "nproc","dt", "tmax"]
data = collect_results(datadir("sims", "dispersion", "frozenfield"), white_list = white_list)

sort!(data, :N)    
transform!(data, :t => ByRow( t -> normalize(fit(Histogram, vec(t), 0:0.001:129), mode = :pdf)) => :tpdf)
transform!(data, :t => ByRow(mean) => :avgt)
transform!(data, :t => ByRow(t->mean(min.(t, 128))) => :avgt128)
transform!(data, :t => ByRow(t->mean(min.(t, 64))) => :avgt64)
transform!(data, :t => ByRow(t->mean(min.(t, 32))) => :avgt32)
transform!(data, :t => ByRow(t->mean(min.(t, 16))) => :avgt16)
transform!(data, :t => ByRow(t->mean(min.(t, 8))) => :avgt8)

ξ=1/3
α =0.5*(2/3-ξ)^(1/2)

# τ pdf and avg
begin
    ξ = 8/9
    α=0
    #equivalent
    #α =0.5*(2/3-ξ)^(1/2)
    
    np = 2
    tmax = 128
    #ncor = 100
    #nindep = 200
    
    df = data[data.ξ .== ξ, :]
    df = df[df.α .== α, :]
    #df = df[df.tmax .== tmax, :]
    #df = df[df.ncor .== ncor, :]
    #df = df[df.nindep .== nindep, :]
    
    plt1 = plot(title = "SDE, ξ = $ξ ,α = $α ", xlabel = "τ^($np)" , ylabel = "PDF_N",xaxis=:linear,yaxis=:linear, xlims = (0,tmax), ylims = (0,2), legend = :topright);
    for i in 1:nrow(df) 
        s= df.nindep[i]*df.nproc[i]
        plot!(Array(range(0,stop=129,length=129000),),  df.tpdf[i].weights, label="sample = $s ,N=$(df.N[i]),dt=$(df.dt[i])") #, c=cm[j][k+1])    
    end
    
    plt2 = plot(title = "", xlabel = "η" , ylabel = "<τ>", xaxis=:log, yaxis=:log, xlims = (1e-4,10), ylims = (1e-2,300), legend = :topright);
    for i in 1:nrow(df)
        plt2 = scatter!([8*pi/df.N[i]], [df.avgt8[i]],markershape=:circle,label="")
        plt2 = scatter!([8*pi/df.N[i]], [df.avgt16[i]],markershape=:x,label="")
        plt2 = scatter!([8*pi/df.N[i]], [df.avgt32[i]],markershape=:square,label="")
        plt2 = scatter!([8*pi/df.N[i]], [df.avgt64[i]],markershape=:diamond,label="")
        plt2 = scatter!([8*pi/df.N[i]], [df.avgt128[i]],markershape=:hexagon,label="")
    end

    layout = @layout [a ; c]
    @show plot(plt1, plt2, layout = (2,1))
end

D0=1
αmean(λ, ξ, α) = (λ^2)/(2 * ((2/(2-(ξ+4*α^2)))*(1-(ξ+4*α^2)^2)/(1+(ξ+4*α^2))) * (D0*(1-((ξ+4*α^2)/2))^2))

αmean(1,8/9,0)