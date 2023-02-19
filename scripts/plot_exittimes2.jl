using DrWatson
@quickactivate("IntermittentFields");
using JLD2
using Statistics
using StatsPlots
using Plots;plotly()
#using StatsBase
#using LinearAlgebra

# theoretical average
D0 = 1
gaussianmean(λ, ξ) = (λ^2)/(2 * ((2/(2-ξ))*(1-ξ^2)/(1+ξ)) * (D0*(1-(ξ/2))^2))

# <τ(λ,ξ)>vsλ  
begin
    plt = plot(title = "<τ(λ,ξ)>", xlabel = "λ", legend = :topleft);
    plot!(xaxis = :log, yaxis = :log);
    plot!(xlims = (1e-2, 2), ylims = (1e-3, 10));
    #plot!(xaxis = :linear, yaxis = :linear);
    #plot!(xlims = (0, 1.2), ylims = (0, 3));
    
    #files = readdir(datadir("sims", "2particles", "non_rescaled_kernel"), join = true)
    files = readdir(datadir("sims", "2particles", "piecewisekernel"), join = true)
    #files = filter(x->occursin("N=8192", x), files)
    files = filter(x->occursin("N=4096", x), files)
    files = filter(x->occursin("dt=0.000156", x), files)

    λs = [1/2^i for i in 0:6]
    cm = cgrad(:jet, length(files), categorical = true)
    
    for (i,f) in enumerate(files)
        g = jldopen(f, "r") 
        t = g["t"]
        t2 = g["t2"]
        t4 = g["t4"]
        t8 = g["t8"]
        t16 = g["t16"]
        t32 = g["t32"]
        t64 = g["t64"]
        
        #data points
        ξ = parse(Float64, split(split(f, "ξ=")[2], ".jld2")[1])
        α = parse(Float64, split(split(f, "α=")[2], "_")[1])
        #if α == 0 
        scatter!(λs, [mean(t), mean(t2), mean(t4), mean(t8), mean(t16), mean(t32), mean(t64)], label = "ξ = $ξ", ms = 4, mc = cm[i])
        
        # theoretical 
        λ = 0:0.001:1
        plot!(λ, gaussianmean.(λ, ξ), label = "ξ = $ξ", lw = 2, c = cm[i])
    end
    display(plt)
end


# <τ(λ,ξ)> vs ξ  
begin
    plt = plot(title = "<τ(λ,ξ)>", xlabel = "ξ", legend = :topleft);
    #plot!(xaxis = :log, yaxis = :log);
    #plot!(xlims = (1e-1, 2), ylims = (1e-3, 10));
    plot!(xaxis = :linear, yaxis = :linear);
    plot!(xlims = (0, 1.2), ylims = (0, 5));
    
    #files = readdir(datadir("sims", "2particles", "non_rescaled_kernel"), join = true)
    files = readdir(datadir("sims", "2particles", "piecewisekernel"), join = true)
    #files = filter(x->occursin("N=8192", x), files)
    files = filter(x->occursin("N=4096", x), files)
    files = filter(x->occursin("dt=0.000156", x), files)

    λs = [1/2^i for i in 0:6]
    cm = cgrad(:rainbow, 7, categorical = true)
    
    for (i,f) in enumerate(files)
        g = jldopen(f, "r") 
        t = g["t"]
        t2 = g["t2"]
        t4 = g["t4"]
        t8 = g["t8"]
        t16 = g["t16"]
        t32 = g["t32"]
        t64 = g["t64"]
        
        #data points
        ξ = parse(Float64, split(split(f, "ξ=")[2], ".jld2")[1])
        #α = parse(Float64, split(split(f, "α=")[2], "_")[1])
        #if α == 0 
        scatter!([ξ], [mean(t)], label="", ms = 4, mc = cm[1])
        scatter!([ξ], [mean(t2)],label="",ms = 4, mc = cm[2])
        scatter!([ξ], [mean(t4)],label="",ms = 4, mc = cm[3])
        scatter!([ξ], [mean(t8)],label="",ms = 4, mc = cm[4])
        scatter!([ξ], [mean(t16)],label="",ms = 4, mc = cm[5])
        scatter!([ξ], [mean(t32)],label="",ms = 4, mc = cm[6])
        scatter!([ξ], [mean(t64)],label="",ms = 4, mc = cm[7])        
        
    end
    # theoretical 
    ξs = 0:0.001:1
    plot!(ξs, gaussianmean.(1, ξs), label = "λ = 1", lw = 2, c = cm[1])
    plot!(ξs, gaussianmean.(1/2, ξs), label = "λ = 1/2", lw = 2, c = cm[2])
    plot!(ξs, gaussianmean.(1/4, ξs), label = "λ = 1/4", lw = 2, c = cm[3])
    plot!(ξs, gaussianmean.(1/8, ξs), label = "λ = 1/8", lw = 2, c = cm[4])
    plot!(ξs, gaussianmean.(1/16, ξs), label = "λ = 1/16", lw = 2, c = cm[5])
    plot!(ξs, gaussianmean.(1/32, ξs), label = "λ = 1/32", lw = 2, c = cm[6])
    plot!(ξs, gaussianmean.(1/64, ξs), label = "λ = 1/64", lw = 2, c = cm[7])
   
    display(plt)
end

