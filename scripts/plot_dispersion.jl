using DrWatson
@quickactivate("IntermittentFields");
using JLD2
using Statistics
using StatsPlots;plotly()
using StatsBase
using LinearAlgebra


# PDF of shape  =========================================================================
#path = "old"
#prefix = "ET3"
path  = "piecewisekernel"
prefix = "ET3l2normR12R13"

begin
    α = 0.0
    plt = plot(title= "pdf of shapes α=$α, Np=4e5",ylims=(0,12),xlims=(0,1), legend = :topleft);
    plot!(yaxis = "PDF_N")
    #plot!(yscale = :log, ylims = (1e-3,10))
    for (j, α) in enumerate([α])
        for (k, i) in enumerate(7:14)
            N = 2^i
            allf = filter(x->occursin(prefix, x), readdir(datadir("sims","zeromodes",path, "dispersion"), join = true))
            allf = filter(x->occursin("N=$N", x), allf)
            allf = filter(x->occursin("α=$α", x), allf)
            
            f1 = jldopen(allf[1], "r")
            #T = f1["t"]
            p = f1["p"]
            #choosing observable
            r1 = vec(abs.(p[1,:,:]-p[2,:,:]))
            r2 = vec(abs.(p[2,:,:]-p[3,:,:]))
            r3 = vec(abs.(p[1,:,:]-p[3,:,:]))
            close(f1)
            if length(allf)>1 
                for i in 2:length(allf)
                    f = jldopen(allf[i], "r")
                    #T = [T f["t"]]
                    p = f["p"]
                    #choosing observables
                    r1 = [r1; vec(abs.(p[1,:,:]-p[2,:,:]))]
                    r2 = [r2; vec(abs.(p[2,:,:]-p[3,:,:]))]
                    r3 = [r3; vec(abs.(p[1,:,:]-p[3,:,:]))] 
                end
            end
            l = length(r1)
            pdf = normalize(fit(Histogram,r1, 0:0.001:1.0), mode=:pdf)
            plot!(Array(range(0,stop=1,length=1000),),  pdf.weights, label="N=$N, correlated") #, c=cm[j][k+1])
            #plot!(title = "Np = $l, α =$α")
        end
        
    end
    plt = @show plt
end

# pdf of τ^(3)  =========================================================================
begin
    α = 0.0
    plt = plot(title= "pdf of τ_3 α=$α, Np=4e5", ylabel = "PDF", xlabel = "τ^(3)")
    #plot!(ylims=(1e-3,10),xlims=(0,10), xscale = :linear, yscale = :log);
    plot!(ylims=(0,4),xlims=(0,10),  yscale = :linear);
    for (j, α) in enumerate([α])
        for (k, i) in enumerate(7:14)
            N = 2^i
            allf = filter(x->occursin("ET3l2normR12R13", x), readdir(datadir("sims","zeromodes","piecewisekernel","dispersion"), join = true))
            allf = filter(x->occursin("N=$N", x), allf)
            allf = filter(x->occursin("α=$α", x), allf)
            
            f1 = jldopen(allf[1], "r")
            t = f1["t"]
            p = f1["p"]
            #choosing observable
            r1 = vec(abs.(p[1,:,:]-p[2,:,:]))
            r2 = vec(abs.(p[2,:,:]-p[3,:,:]))
            r3 = vec(abs.(p[1,:,:]-p[3,:,:]))
            t =  vec(t)
            close(f1)
            if length(allf)>1 
                for i in 2:length(allf)
                    f = jldopen(allf[i], "r")
                    t = [t; vec(f["t"])]
                    p = f["p"]
                    #choosing observables
                    r1 = [r1; vec(abs.(p[:,1,:]-p[:,2,:]))]
                    r2 = [r2; vec(abs.(p[:,2,:]-p[:,3,:]))]
                    r3 = [r3; vec(abs.(p[:,1,:]-p[:,3,:]))] 
                end
            end
            
            pdf = normalize(fit(Histogram, t, 0:0.1:10), mode=:pdf)
            plot!(Array(range(0,stop=10,length=100),),  pdf.weights, label="N=$N") #, c=cm[j][k+1])
        end
    end
    plt = @show plt
end
