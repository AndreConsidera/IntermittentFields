using Distributed
nproc = 8
procs = addprocs(nproc)
@everywhere using DrWatson
@everywhere @quickactivate "IntermittentFields";
@everywhere include(srcdir("IntermittentFields.jl"))
@everywhere using .IntermittentFields
@everywhere using FFTW

using LinearAlgebra
using Statistics
@everywhere using Distributions
@everywhere using Interpolations
using Plots;plotly(size=(800,650))
using DataFrames
using Statistics
using StatsBase


@everywhere function quenchedgmc(nindep::Integer,dt::Real, N::Real ,α::Real, ξ::Real, tmax::Real, beta::Real, FFTW_effort::UInt32, c::Function, np::Integer, ncor::Integer, sizefunc::Function, λ::Real)
    dx = 2 * pi/N    
    eta = 2 * dx
    r0 = 2 * eta
    r = [-pi + i * dx for i in 0:N-1]
    kappa = eta^beta
    
    #instantiating kernels and FFT
    ker = CovarianceKernel(r, eta, CovarianceCorrelation(c), ξ, false)
    logker = SingularKernel(r, eta)
    Pinv = plan_irfft(copy(ker.Lk), length(ker.r); flags = FFTW_effort, timelimit = Inf)
    Pforward = plan_rfft(copy(ker.r); flags = FFTW_effort, timelimit = Inf);
    
    push!(r, pi)
    maxiter = ceil(tmax/dt);
    exittime = zeros(ncor, nindep);
    finalp = zeros(np, ncor, nindep);
    for k in 1:nindep
        #initial particles
        p =  zeros(np, ncor)
        p[1,:] = rand(Uniform(-pi, pi), ncor)
        if np > 1
            for l in 2:np  
                p[l , :] = p[l - 1,:] .+ r0
            end
        end
        t = zeros(ncor)
        
        # fix the gmc during time evolution 
        g2 = UnitaryWhiteNoise(div(length(ker.r), 2) + 1)
        for i in 1:maxiter
            R = sizefunc(p)
            if minimum(R) >= λ break end
            t[R .< λ] = t[R .< λ] .+ dt

            g1 = UnitaryWhiteNoise(div(length(ker.r), 2) + 1)
            # here we use Pforward and Pinv
            noise = GmcNoise(logker, g1, g2, α, Pforward, Pinv)
            u = realization(ker, noise, Pinv)
            uinterp = linear_interpolation(r, push!(u, u[1]), extrapolation_bc = Periodic())
            p[:, R .< λ] .= p[:, R .< λ] .+ uinterp.(p[:, R .< λ]) .* dt.^0.5 .+ (2 .* kappa .* dt).^0.5 .* rand(Normal(0,1), np, size(p[:, R .< λ], 2) )
        end
        exittime[:, k] = t
        finalp[:, :, k] = p
    end
    return exittime, finalp
end


# TRAJECTORY
tmax=10
Np=1
dt=tmax/2^13
maxiter=Int(ceil(tmax/dt))
t= zeros(maxiter)
x=zeros(maxiter, 4)

N=2^14
dx = 2 * pi/N 
r = [-pi + i * dx for i in 0:N-1]
eta = 4*pi/N;
ξ = 1/3
logker = SingularKernel(r, eta) 
σ_sq = 2 * sum(abs.(logker.Lk[2:end]).^2)

#g2 = UnitaryWhiteNoise(div(length(logker.r), 2) + 1)

α = 0.0
gmc = exp.(α .* realization(logker, g2) .- α^2 .* σ_sq)
p3=plot(r,gmc,xaxis=[-pi,pi],label="α=$α",yaxis=:linear)
pwker = CovarianceKernel(r, eta, CovarianceCorrelation(piecewisekernel), ξ, false);
push!(r, pi)

for k in 1:Np  
    x0=rand(Uniform(0,1),1)[1]
    #y0=x0+eta
    x[1,1] = 0
    x[1,2] = 0 + eta
    x[1,3] = -2
    x[1,4] = -2 +eta
    kappa=eta
    for i in 1:maxiter-1
        t[i+1]= t[i] + dt
        g1 = UnitaryWhiteNoise(div(length(pwker.r), 2) + 1)
        noise = GmcNoise(logker, g1, g2, α)
        
        u = realization(pwker, noise)
        uinterp = linear_interpolation(r, push!(u, u[1]), extrapolation_bc = Periodic())
        x[i + 1, :] .= x[i, :] .+ uinterp.(x[i,:]) .* dt.^0.5 .+ (2 .* kappa .* dt).^0.5 .* rand(Normal(0,1), size(x, 2))
    end
end

p3=plot!(x,t,xaxis=[-pi,pi])

p2=plot!(x,t,xaxis=[-pi,pi])

pl=plot!(y,t,xaxis=[-pi,pi])

layout = @layout [a;c ]
plot(p1,p2,layout=layout)

# EXIT TIMES

# sizefunc
@everywhere R12R13(p::AbstractArray{<:Real}) = @. sqrt((p[1,:]-p[2,:])^2 + (p[1,:]-p[3,:])^2)
@everywhere R12(p::AbstractArray{<:Real}) = @. sqrt((p[1,:]-p[2,:])^2)

# all parameters
d = Dict(
"np" => [2], 
"ncor" => [1],
"nindep" => [100],
"α" => [0.6],
"ξ" => [2/3],
"kernel" => piecewisekernel,
"tmax" => 10,
"λ" => [1],
"N" => 2^10,
"dt"=> 1/6400,
"beta" => 1.0,
"planning_effort" => FFTW.MEASURE,
"sizefunc" => [@onlyif("np" == 2, R12), @onlyif("np" == 3, R12R13)], 
)

dicts = dict_list(d);
    
function makesim(d::Dict)
    @unpack np, ncor, nindep, α, ξ, kernel, tmax, λ, N, dt, beta, planning_effort, sizefunc = d
    t = zeros(ncor, nindep, nproc);
    p = zeros(np, ncor, nindep, nproc);
    a = Array{Future}(undef, nproc)
    for i in 1:nproc
        a[i] = @spawnat procs[i]  quenchedgmc(nindep, dt, N, α, ξ, tmax, beta, planning_effort, kernel, np, ncor, sizefunc, λ)
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
    @tagsave(datadir("sims", "dispersion", "quenchedgmc", savename("quenched", d, "jld2")), fulldict, safe = DrWatson.readenv("DRWATSON_SAFESAVE", true))
end        



# data analysis
using DrWatson
@quickactivate "IntermittentFields";
using LinearAlgebra
using Statistics
using Plots;plotly(size=(800,650))
using DataFrames
using Statistics
using StatsBase
white_list = ["N", "α", "p", "t", "ξ", "ncor", "nindep","np", "nproc", "λ","dt", "tmax"]
data = collect_results(datadir("sims", "dispersion", "quenchedgmc"), white_list = white_list)

sort!(data, :N)    
transform!(data, :t => ByRow( t -> normalize(fit(Histogram, vec(t), 0:0.001:65), mode = :pdf)) => :tpdf)
transform!(data, :t => ByRow(mean) => :avgt)
transform!(data, :t => ByRow(t->mean(min.(t, 64))) => :avgt64)
transform!(data, :t => ByRow(t->mean(min.(t, 32))) => :avgt32)
transform!(data, :t => ByRow(t->mean(min.(t, 16))) => :avgt16)
transform!(data, :t => ByRow(t->mean(min.(t, 8))) => :avgt8)


# τ pdf and avg
begin
    ξ = 1/9
    α = 0.3
    np = 2
    tmax = 64
    #ncor = 100
    #nindep = 200
    
    df = data[data.ξ .== ξ, :]
    df = df[df.α .== α, :]
    #df = df[df.tmax .== tmax, :]
    #df = df[df.ncor .== ncor, :]
    #df = df[df.nindep .== nindep, :]
    
    plt1 = plot(title = "fields,ξ = $ξ ,α = $α ", xlabel = "τ^($np)" , ylabel = "PDF_N",xaxis=:linear,yaxis=:linear, xlims = (0,tmax), ylims = (0,2), legend = :topright);
    for i in 1:nrow(df) 
        s= df.ncor[i]*df.nindep[i]*df.nproc[i]
        plot!(Array(range(0,stop=65,length=65000),),  df.tpdf[i].weights, label="sample = $s ,N=$(df.N[i]), λ=$(df.λ[i]),dt=$(df.dt[i])") #, c=cm[j][k+1])    
    end

    plt2 = plot(title = "", xlabel = "η" , ylabel = "<τ^($np)>", xaxis=:log, yaxis=:log, xlims = (1e-4,10), ylims = (1e-2,300), legend = :topright);
    for i in 1:nrow(df)
        plt2 = scatter!([8*pi/df.N[i]], [df.avgt8[i]],markershape=:circle,label="")
        plt2 = scatter!([8*pi/df.N[i]], [df.avgt16[i]],markershape=:x,label="")
        plt2 = scatter!([8*pi/df.N[i]], [df.avgt32[i]],markershape=:square,label="")
        plt2 = scatter!([8*pi/df.N[i]], [df.avgt64[i]],markershape=:diamond,label="")
    end
    layout = @layout [a ; c]
    @show plot(plt1, plt2, layout = layout)
end


# Stochasticity Phase Diagram
# (ξ, α)
stoch = [
    (2/3, 0.0),
    (2/3, 0.1),
    (1/3, 0.0),
    (1/3, 0.1),
    (1/3, 0.2),
    (4/9, 0.0),
    (4/9, 0.1),
    (4/9, 0.2),
    (1/9, 0.0),
    (2/9, 0.0),
    (1/9, 0.1),
    (2/9, 0.1),
    (1/9, 0.2),
    (2/9, 0.2),
    (2/3, 0.2),
    (5/9, 0.0),
    (5/9, 0.1),
    (5/9, 0.2),
]
nostoch = [
    (2/3, 0.4),
    (2/3, 0.5),
    (2/3, 0.6),
    (1/3, 0.6),
    (1/3, 0.5),
    (1/3, 0.4),
    (4/9, 0.4),
    (4/9, 0.5),
    (4/9, 0.6),
    (1/9, 0.3),
    (2/9, 0.3),
    (1/9, 0.4),
    (1/9, 0.5),
    (1/9, 0.6),
    (2/9, 0.4),
    (2/9, 0.5),
    (2/9, 0.6),
]
plt=scatter(title = "Stochasticity Phase Diagram", xlims = [0, 1], ylims = [-0.01, 0.7], xlabel = "ξ", ylabel = "α");
for x in stoch 
    scatter!([x[1]], [x[2]], markershape=:square, markercolor = :black, ms=10, label = "")
end
for x in nostoch 
    scatter!([x[1]], [x[2]], markershape=:square, markercolor = :white, ms=10, label = "")
end
@show plt