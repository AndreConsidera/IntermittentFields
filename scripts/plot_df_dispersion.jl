using DrWatson
@quickactivate("IntermittentFields");
using DataFrames
using Statistics
using StatsPlots;plotly()
using StatsBase
using LinearAlgebra

white_list = ["N", "α", "p", "t"]
data = collect_results(datadir("sims","zeromodes","piecewisekernel", "dispersion"), white_list = white_list)
sort!(data, :N)
transform!(data, :p => ByRow(p -> vec(abs.(p[1,:,:]-p[2,:,:]))) => :r12)
transform!(data, :r12 => ByRow( r12 -> normalize(fit(Histogram, r12, 0:0.001:1.0), mode = :pdf)) => :r12pdf)
transform!(data, :t => ByRow( x -> normalize(fit(Histogram, vec(x), 0:0.1:10), mode = :pdf)) => :tpdf)

# pdf of shape
begin 
    α = 0.0
    df = data[data.α .== α,:]    
    
    plt = plot(title = "dataframe", ylabel = "PDF_N", xlims = (0,1), ylims = (0,12), legend = :topleft);
    for i in 1:nrow(df) 
        plot!(Array(range(0,stop=1,length=1000),),  df.r12pdf[i].weights, label="N=$(df.N[i]), α = $(df.α[i]), correlated") #, c=cm[j][k+1])
    end
    @show plt
end

# τ pdf
begin
    α = 0.0
    df = data[data.α .== α,:]
    
    plt = plot(title = "dataframe", xlabel = "τ^(3)" , ylabel = "PDF_N", xlims = (0,10), ylims = (0,4), legend = :topright);
    for i in 1:nrow(df) 
        plot!(Array(range(0,stop=10,length=100),),  df.tpdf[i].weights, label="N=$(df.N[i]),$(df.α[i]), correlated") #, c=cm[j][k+1])
    end
    @show plt
end
