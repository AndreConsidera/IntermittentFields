using DrWatson
@quickactivate("IntermittentFields");
using JLD2
using Statistics
using StatsPlots
#using StatsBase
#using LinearAlgebra

key = "dt=0.000156"
key = "dt=0.000625"
key = "dt=0.0025"
key = "dt=0.01"

allf = filter(x->occursin(key,x), readdir(datadir("sims", "old")))
allf = filter(x->occursin("rescaled", x), allf)
allf = filter(x->occursin("beta=1.0", x), allf)
allf = filter(x->occursin("ξ=0.6", x), allf)
allf = filter(x->occursin("α=0.25", x), allf)

f1 = jldopen(datadir("sims","old",allf[1]),"r")
T = f1["ET"]
fclose(f1)

for i in 2:length(allf)
    f = jldopen(datadir("sims", "old",allf[i]),"r")
    T = [T f["ET"]]
end

mean(T)
var(T)
density(reshape(T,length(T)); bandwidth = 1e-3, trim = false)

# theoretical average
h=0.3
D0=1
xi=2*h
D=D0*(1-h)^2
d_e= (2/(2-xi))*(1-xi^2)/(1+xi)
ave=1/(2*d_e*D)
nu=d_e/2-1
va=1/(16*D^2*(nu+2)*(nu+1)^2)

#bins = collect(0:0.01:10)
#h = fit(Histogram, reshape(T,length(T)), bins; closed = :right)
#plot(h)
