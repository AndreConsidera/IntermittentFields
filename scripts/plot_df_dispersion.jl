using DrWatson
@quickactivate("IntermittentFields");
include(srcdir("IntermittentFields.jl"))
using .IntermittentFields
using DataFrames
using Statistics
using StatsPlots;plotly()
using StatsBase
using LinearAlgebra

white_list = ["N", "α", "p", "t", "ncor", "nindep", "ξ"]
data = collect_results(datadir("sims","zeromodes","piecewisekernel", "dispersion"), white_list = white_list)

function f1(θ) 
    θ[θ .< pi/2]=θ[θ .< pi/2] .+ pi
    return θ
end

function f2(θ) 
    θ[θ .> pi/2] .= θ[θ .> pi/2] .- pi 
    return θ
end
function f3(θ) 
    θ[(pi/4 .< θ .<= pi/2)] .= pi/2 .- θ[(pi/4 .< θ .<= pi/2)]
    return θ
end
function f4(θ) 
    θ[(-pi/2 .< θ .<= -pi/4)] .= pi/2 .- θ[(-pi/2 .< θ .<= -pi/4)]
    return θ
end


sort!(data, :N)
transform!(data, :p => ByRow(p -> vec(normtorus.(p[1,:,:]-p[2,:,:]))) => :r12)
transform!(data, :p => ByRow(p -> vec(normtorus.(p[1,:,:]-p[3,:,:]))) => :r13)

transform!(data, [:r12, :r13] => ByRow((r12, r13) -> @.angle(complex(r12, r13))) => :θ)

transform!(data, :θ => ByRow(f1)  => :θ)
transform!(data, :θ => ByRow(f2)   => :θ)
transform!(data, :θ => ByRow(f3)  => :θ)
transform!(data, :θ => ByRow(f4)  => :θ)

transform!(data, :θ => ByRow( θ -> normalize(fit(Histogram, θ, -pi/4:0.001:pi/4), mode = :pdf)) => :θpdf)
transform!(data, :r12 => ByRow( r12 -> normalize(fit(Histogram, r12, -0.1:0.001:1.2), mode = :pdf)) => :r12pdf)
transform!(data, :t => ByRow( t -> normalize(fit(Histogram, vec(t), 0:0.01:10), mode = :pdf)) => :tpdf)
αs = sort(unique(data.α))


# pdf of shape
begin 

    α = 0.6
    df = data[data.α .== α,:]    
    
    plt = plot(title = "Correlated α = $α", ylabel = "PDF_N", xlims = (-0.1,1.2), ylims = (0,20), legend = :topleft);
    for i in 1:nrow(df) 
        plot!(Array(range(0,stop=1.1,length=1300),),  df.r12pdf[i].weights, label="N=$(df.N[i]), correlated") #, c=cm[j][k+1])
    end
    @show plt
end

# τ pdf

begin
    α = 0.0
    Nmax = 512
    df = data[data.α .== α, :]
    df =df[df.N .<= Nmax, :]
    df=df[df.nindep.==256,:]
    plt = plot(title = "correlated α = $α", xlabel = "τ^(3)" , ylabel = "PDF_N", xlims = (0,3), ylims = (0,4), legend = :topright);
    for i in 1:nrow(df) 
        plot!(Array(range(0,stop=10,length=1000),),  df.tpdf[i].weights, label="N=$(df.N[i]), ncor=$(df.ncor[i]), nindep=$(df.nindep[i])") #, c=cm[j][k+1])
    end
    @show plt
end

# θpdf
begin
    α = 0.0
    Nmax = 2*4096
    df = data[data.α .== α, :]
    df =df[df.N .== Nmax, :]
    plt = plot(title = "correlated α = $α", xlabel = "τ^(3)" , ylabel = "PDF_N", xlims = (-pi/4,20*pi/4), ylims = (0,20), legend = :topright);
    for i in 1:nrow(df) 
        plot!(Array(range(0,stop=10,length=1570),),  df.θpdf[i].weights, label="N=$(df.N[i]), ncor=$(df.ncor[i]), nindep=$(df.nindep[i])") #, c=cm[j][k+1])
    end
    @show plt
end

