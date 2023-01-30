using DrWatson
@quickactivate("Intermittency Paradox");
using JLD2
using Statistics
using StatsPlots;plotly()
using StatsBase
using LinearAlgebra

dt=0.000156
α=0.0

#reference values
allf = filter(x->occursin("ET3l2normR12R13", x), readdir(datadir("sims","zeromodes")))
allf = filter(x->occursin("N=128", x), allf)
allf = filter(x->occursin("dt=$dt", x), allf)
#allf = filter(x->occursin("beta=1.0", x), allf)
#allf = filter(x->occursin("ξ=0.6", x), allf)
allf = filter(x->occursin("α=$α", x), allf)

f1 = jldopen(datadir("sims", "zeromodes", allf[1]), "r")
T = f1["t"]
p = f1["p"]
#choosing observable
R1 = vec(abs.(p[:,1,:]-p[:,2,:]))
R2 = vec(abs.(p[:,2,:]-p[:,3,:]))
R3 = vec(abs.(p[:,1,:]-p[:,3,:]))
close(f1)
#=    
for i in 2:length(allf)
    f = jldopen(datadir("sims", allf[i]),"r")
    T = [T f["ET"]]
    p = f["p"]
    #choosing observables
    R1 = [R1; vec(abs.(p[:,1,:]-p[:,2,:]))]
    R2 = [R2; vec(abs.(p[:,2,:]-p[:,3,:]))]
    R3 = [R3; vec(abs.(p[:,1,:]-p[:,3,:]))] 
end
=#
pdf_ref = normalize(fit(Histogram, R1, 0:0.001:1), mode=:pdf)
 

#plt = density(yaxis=:log, xaxis=:log, title="min(r1,r2,r3),α=$α, Np=6400,dt=$dt",ylims=(1e-3,20),xlims=(1e-3,1.2));
#plt = density(title="min(r1,r2),α=$α, Np=6400,dt=$dt",ylims=(0,3),xlims=(-1.5,1.5));
plt = plot(title="α=$α, Np=6400,dt=$dt",ylims=(-3,3),xlims=(0,1));

#other values
for i in 8:9
    N = 2^i
    allf = filter(x->occursin("ET3l2normR12R13", x), readdir(datadir("sims", "zeromodes")))
    allf = filter(x->occursin("N=$N", x), allf)
    allf = filter(x->occursin("dt=$dt", x), allf)
    #allf = filter(x->occursin("beta=1.0", x), allf)
    #allf = filter(x->occursin("ξ=0.6", x), allf)
    allf = filter(x->occursin("α=$α", x), allf)
    
    f1 = jldopen(datadir("sims", "zeromodes", allf[1]), "r")
    T = f1["t"]
    p = f1["p"]
    #choosing observable
    r1 = vec(abs.(p[:,1,:]-p[:,2,:]))
    r2 = vec(abs.(p[:,2,:]-p[:,3,:]))
    r3 = vec(abs.(p[:,1,:]-p[:,3,:]))
    close(f1)
    
    if length(allf)>1 
        for i in 2:length(allf)
            f = jldopen(datadir("sims", "zeromodes", allf[i]),"r")
            T = [T f["t"]]
            p = f["p"]
        #choosing observables
            r1 = [r1; vec(abs.(p[:,1,:]-p[:,2,:]))]
            r2 = [r2; vec(abs.(p[:,2,:]-p[:,3,:]))]
            r3 = [r3; vec(abs.(p[:,1,:]-p[:,3,:]))] 
        end
    end
    
    # min(r1,r2,r3)
    #density!(min.(r1, r2) - min.(R1, R2), label="N=$N"; bandwidth = 1e-3, trim = false)
    #density!(r1, label="N=$N"; bandwidth = 1e-3, trim = false)
    
    pdf = normalize(fit(Histogram,r1, 0:0.001:1.0), mode=:pdf)

    plot!(Array(range(0,stop=1,length=1000),), (N/128)^(0.5) *(pdf.weights - pdf_ref.weights),label="N=$N")
end
plt =@show plt

savefig("zeromode_α=$α.png")

# T
mean(T)
var(T)
density(title="Exit time:|r1|+|r2|>1, α=$α, Np=6400,dt=$dt",ylims=(0,1.5),xlims=(0,5));
density!(reshape(T, length(T)); bandwidth = 1e-3, trim = false)

#savefig(plotsdir("ET3_α=$α.png"))

density(title="r1, Np=6400,dt=$dt",ylims=(0,12),xlims=(0,1.2));
density!(reshape(r1, length(r1)),label="α=$α"; bandwidth = 1e-3, trim = false)
#savefig(plotsdir("r1_α=$α.png"))


density(title="r1/(r1+r2), Np=6400,dt=$dt",ylims=(0,17),xlims=(0,1.2));
density!(reshape(r1./(r1+r2), length(r1)),label="α=$α"; bandwidth = 1e-3, trim = false)
savefig(plotsdir("shape_α=$α.png"))




a=vec(p[:,1,:]-p[:,2,:])

ppp =reshape(pp,(3,800))

ppp[1,:]-ppp[2,:]