using DrWatson
@quickactivate("IntermittentFields");
using JLD2
using Statistics
using StatsPlots;plotly()
using StatsBase
using LinearAlgebra


# PDF of shape  =========================================================================
#path = "old"
path  = "piecewisekernel"
prefix = "ET3l2normR12R13"
#prefix = "ET3"

begin
    α = 0.6
    plt = plot(title= "pdf of shapes α=$α, Np=4e5",ylims=(0,12),xlims=(0,1));
    plot!(yaxis = "PDF_N")
    #plot!(yscale = :log, ylims = (1e-3,10))
    for (j, α) in enumerate([α])
        for (k, i) in enumerate(9:11)
            N = 2^i
            allf = filter(x->occursin(prefix, x), readdir(datadir("sims","zeromodes",path, "dispersion"), join = true))
            allf = filter(x->occursin("N=$N", x), allf)
            allf = filter(x->occursin("α=$α", x), allf)
            
            f1 = jldopen(datadir("sims", "zeromodes", path, "dispersion", allf[1]), "r")
            #T = f1["t"]
            p = f1["p"]
            #choosing observable
            r1 = vec(abs.(p[1,:,:]-p[2,:,:]))
            r2 = vec(abs.(p[2,:,:]-p[3,:,:]))
            r3 = vec(abs.(p[1,:,:]-p[3,:,:]))
            close(f1)
            if length(allf)>1 
                for i in 2:length(allf)
                    f = jldopen(datadir("sims", "zeromodes", path ,"dispersion", allf[i]),"r")
                    #T = [T f["t"]]
                    p = f["p"]
                    #choosing observables
                    r1 = [r1; vec(abs.(p[1,:,:]-p[2,:,:]))]
                    r2 = [r2; vec(abs.(p[2,:,:]-p[3,:,:]))]
                    r3 = [r3; vec(abs.(p[1,:,:]-p[3,:,:]))] 
                end
            end
            
            pdf = normalize(fit(Histogram,r1, 0:0.001:1.0), mode=:pdf)
            plot!(Array(range(0,stop=1,length=1000),),  pdf.weights, label="N=$N, correlated") #, c=cm[j][k+1])
        end
        
    end
    plt = @show plt
end
