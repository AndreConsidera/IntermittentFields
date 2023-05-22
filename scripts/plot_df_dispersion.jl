using DrWatson
@quickactivate("IntermittentFields");
include(srcdir("IntermittentFields.jl"))
using .IntermittentFields
using DataFrames
using Statistics
using StatsPlots;plotly(size=(800,650))
using StatsBase
using LinearAlgebra
using PyFormattedStrings	

white_list = ["N", "α", "p", "t", "ξ", "ncor", "nindep","np", "nproc", "λ","dt"]
data = collect_results(datadir("sims", "dispersion","data_for_paper"), white_list = white_list , rinclude = [r"ET3l2normR12R13"])

function f(θ) 
    θ[θ .< pi/2]=θ[θ .< pi/2] .+ pi
    θ[θ .> pi/2] .= θ[θ .> pi/2] .- pi 
    θ[(pi/4 .< θ .<= pi/2)] .= pi/2 .- θ[(pi/4 .< θ .<= pi/2)]
    θ[(-pi/2 .< θ .<= -pi/4)] .= pi/2 .- θ[(-pi/2 .< θ .<= -pi/4)]
    return θ
end

sort!(data, :N)
transform!(data, :p => ByRow(p -> vec(normtorus.(p[1,:,:,:]-p[2,:,:,:]))) => :r12)
transform!(data, :p => ByRow(p -> vec(normtorus.(p[1,:,:,:]-p[3,:,:,:]))) => :r13)
transform!(data, [:r12, :r13] => ByRow((r12, r13) -> @.angle(complex(r12, r13))) => :θ)
transform!(data, :θ => ByRow(f)  => :θ)
transform!(data, :θ => ByRow( θ -> normalize(fit(Histogram, θ, -pi/4:0.001:pi/4), mode = :pdf)) => :θpdf)
transform!(data, :r12 => ByRow( r12 -> normalize(fit(Histogram, r12, -0.1:0.001:1.2), mode = :pdf)) => :r12pdf)
transform!(data, :t => ByRow( t -> normalize(fit(Histogram, vec(t), 0:0.001:10), mode = :pdf)) => :tpdf)
transform!(data, :t => ByRow(mean) => :avgt)
αs = sort(unique(data.α))


# pdf of shape
begin 
    α = 0.0
    df = data[data.α .== α,:] 
    plt = plot(title = "Correlated α = $α", ylabel = "PDF_N", xlims = (0,1), ylims = (0,14), legend = :topleft);
    for i in 1:nrow(df) 
        plot!(Array(range(0,stop=1.1,length=1300),),  df.r12pdf[i].weights, label="N=$(df.N[i]), correlated") #, c=cm[j][k+1])
    end
    @show plt
end

# τ pdf and avg
begin
    α = 0.0
    ξ = 1/3
    #dt= 1e-4
    N = 2^10
    λ = 1
    np = 2
    ncor = 100
    df = data[data.np .== np, :]
    df = df[df.λ .== λ, :]
    df = df[df.ξ .== ξ, :]
    df = df[df.ncor .== ncor, :]
    #df = df[df.N .== N, :]
    D0 = 1
    gaussianmean(λ, ξ) = (λ^2)/(2 * ((2/(2-ξ))*(1-ξ^2)/(1+ξ)) * (D0*(1-(ξ/2))^2))
    
    plt1 = plot(title = "", xlabel = "τ^($np)" , ylabel = "PDF_N",xaxis=:linear,yaxis=:log, xlims = (0,10), ylims = (0,2), legend = :topright);
    for i in 1:nrow(df) 
        s= df.ncor[i]*df.nindep[i]*df.nproc[i]
        plot!(Array(range(0,stop=10,length=10000),),  df.tpdf[i].weights, label="sample = $s ,N=$(df.N[i]), λ=$(df.λ[i]),dt=$(df.dt[i])") #, c=cm[j][k+1])    
    end
    
    plt2 = plot(title = "", xlabel = "η" , ylabel = "<τ^($np)>", xaxis=:log, yaxis=:linear, xlims = (1e-4,10), ylims = (0,3), legend = :topright);
    for i in 1:nrow(df)
        plt2 = scatter!([8*pi/df.N[i]], [df.avgt[i]])
    end
    layout = @layout [a ; c]
    @show plot(plt1, plt2, layout = layout)
end

# θpdf
begin
    α = 0.0
    Nmax = 2*4096
    df = data[data.α .== α, :]
    df =df[df.N .== Nmax, :]
    plt = plot(title = f"correlated α = {α:0.2f}, ξ = {ξ:0.2f}", xlabel = "τ^(3)" , ylabel = "PDF_N", xlims = (-pi/4,20*pi/4), ylims = (0,20), legend = :topright);
    for i in 1:nrow(df) 
        plot!(Array(range(0,stop=10,length=1570),),  df.θpdf[i].weights, label="N=$(df.N[i]), ncor=$(df.ncor[i]), nindep=$(df.nindep[i])") #, c=cm[j][k+1])
    end
    @show plt
end


# <τ> x η
begin
    α = 0.0
    ξ = 1/3
    dt= 1e-5 #1/(4*16*100)
    #N = 2^18
    df = data[data.np .== 2, :]
    df = df[df.ξ .== ξ, :]
    df = df[df.dt .== dt, :]
    df = df[df.N .== N, :]
    df = df[df.α .== α, :]
    D0 = 1
    gaussianmean(λ, ξ) = (λ^2)/(2 * ((2/(2-ξ))*(1-ξ^2)/(1+ξ)) * (D0*(1-(ξ/2))^2))
    
    plt = plot(title = f" α = {α:0.2f}, ξ = {ξ:0.3f}" , ylabel = "<τ(2)>",xlabel="η", xlims = (0,1.2), ylims = (2.5,5.5), legend = :topleft);
    #plot!(inset = (1, bbox(0.1,0.1, 0.5,0.5)))
    #plot!(subplot = 2 ,yscale=:log,xscale=:log, xlims = (1e-3,5.2),ylim=(1e-6,10))
    #plot!([i for i in 0:0.0001:1],gaussianmean.([i for i in 0:0.0001:1], ξ), subplot=2, label = "")
   #plot!([i for i in 0:0.00001:1],gaussianmean.([i for i in 0:0.00001:1], ξ) ./([i for i in 0:0.00001:1]).^2);
    
    for i in 1:nrow(df) 
        scatter!([1/df.N[i]],  [df.avgt[i]], label="") #, c=cm[j][k+1])
        #scatter!([df.λ[i]],  [df.avgt[i]]/(df.λ[i])^2, label= "ncor=$(df.ncor[i]),N = $(df.N[i]),dt = $(df.dt[i])") #, c=cm[j][k+1])
    end
    @show plt

end

