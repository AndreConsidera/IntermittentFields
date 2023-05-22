using Distributed
nproc = 8
procs = addprocs(nproc)
@everywhere using DrWatson
@everywhere @quickactivate  
@everywhere include(srcdir("IntermittentFields.jl"))
@everywhere using .IntermittentFields
@everywhere using LinearAlgebra
using Plots;plotlyjs(size=(800,650))
using DataFrames
using Statistics
using StatsBase
@everywhere using Distributions
@everywhere using Interpolations
using Random
@everywhere using DifferentialEquations
using DifferentialEquations.EnsembleAnalysis


@everywhere function sdesolver(nindep::Integer,dt::Real, N::Real ,α::Real, ξ::Real, tmax::Real, beta::Real)

    dx = 2 * pi/N    
    eta = 2 * dx
    kappa =eta^beta
    r = [-pi + i * dx for i in 0:N-1]
    logker = SingularKernel(r, eta)
    push!(r, pi)
    σ_sq = 2 * sum(abs.(logker.Lk[2:end]).^2)

    
    #setting environment
    g2 = UnitaryWhiteNoise(div(length(logker.r), 2) + 1)
    gmc = exp.(α .* realization(logker, g2) .- α^2 .* σ_sq)
    gmcinterp = linear_interpolation(r, push!(gmc, gmc[1]), extrapolation_bc = Periodic())    
    p=[gmcinterp]
    prob_func = let p = p
        (prob, i, repeat) -> begin
            g2 = UnitaryWhiteNoise(div(length(logker.r), 2) + 1)
            gmc = exp.(α .* realization(logker, g2) .- α^2 .* σ_sq)
            gmcinterp = linear_interpolation(r, push!(gmc, gmc[1]), extrapolation_bc = Periodic())
        
            remake(prob, p = [gmcinterp])
        end
    end
    
    function piecewisekernel(x, y, η, ξ)
        r=x-y
        @.((1. - normevans(r,η)^ξ) * (r<=1)) 
    end

    function f(du,u, p, t) 
        du[1]=0.
        du[2]=0.
    end

    function g(du, u, p, t)
        L_sq = [piecewisekernel(x, y, eta, ξ) for x in u, y in u]
        if isposdef(L_sq)
            L, U = cholesky(L_sq)
        end
        gmcinterp = p[1]
        du[1,1] = L[1,1]*gmcinterp(u[1])
        du[1,2] = L[1,2]*gmcinterp(u[2])
        du[2,1] = L[2,1]*gmcinterp(u[1])
        du[2,2] = L[2,2]*gmcinterp(u[2])
        #kappa
        du[1,3] = (2*kappa)^0.5
        du[2,3] = (2*kappa)^0.5
    end

    function terminate_condition(u, t, integrator) 
        abs(u[1]-u[2]) - 1
    end
    
    function eta_boundary_condition(u, t, integrator) 
        abs(u[1]-u[2]) < eta
    end
    
    terminate_affect!(integrator) = terminate!(integrator)
    function eta_boundary_affect!(integrator)  
        integrator.u[1]=0.5*(integrator.u[1]+integrator.u[2])+eta/2
        integrator.u[2]=0.5*(integrator.u[1]+integrator.u[2])-eta/2
    end
    
    cb1 = ContinuousCallback(terminate_condition, terminate_affect!)
    cb2 = SciMLBase.DiscreteCallback(eta_boundary_condition, eta_boundary_affect!)
    cbs=CallbackSet(cb1, cb2)

    u0 = [0.,eta]
    tspan = (0.0, tmax)
    prob = SDEProblem(f, g, u0, tspan, p, noise_rate_prototype = zeros(2, 3))
    ensembleprob = EnsembleProblem(prob, output_func = (sol, i) -> (sol.t[end], false), prob_func = prob_func)
    global sol = solve(ensembleprob, EM(), EnsembleDistributed(); dt = dt, callback = cbs, trajectories = nindep)
    return sol
end

# all parameters
d = Dict{String, Any}( 
"nindep" => 100,
"α" => 0.4,
"ξ" => 2/3,
"tmax" => 128,
"N" => 2^14,
"dt"=> 1e-4,
"beta" => 1.0,
"placeholder"=>"togetcorrecttype "
)

dicts = dict_list(d);
    
function makesim(d::Dict)
    @unpack nindep, α, ξ, tmax, N, dt, beta = d
    sol=sdesolver(nindep::Integer,dt::Real, N::Real ,α::Real, ξ::Real, tmax::Real, beta::Real)
    
    fulld = copy(d)
    fulld["sol"] = sol
    fulld["nproc"] = nproc
    return fulld
end


for (i, d) in enumerate(dicts) 
    println("sim number=$i/$(length(dicts))")
    println("parameters:", savename(d),"_nproc=$nproc")
    global fulldict = makesim(d)
   # @tagsave(datadir("sims", "dispersion", "sdesolver", savename("solver", d, "jld2")), fulldict, safe = DrWatson.readenv("DRWATSON_SAFESAVE", true))
end     

mean(fulldict["sol"].u)
@show sol.u
#=
summ = EnsembleSummary(sol)
plot(summ,error_style = :none)
plot(summ; idxs = 1,ci_type=:SEM)
plot!(sol[1].t,1.96sol[1].t.^(1/2))
size(sol)

ts=sol.u[1].u
plot(r,gmcinterp(r))
plot(sol.u[1])
m_series, v_series = timeseries_steps_meanvar(sol)
plot(v_series)

#sol = solve(prob, SRIW1(), dt = dt, adaptive = false)

D0 = 1
αmean(λ, ξ, α) = (λ^2)/(2 * ((2/(2-(ξ+4*α^2)))*(1-(ξ+4*α^2)^2)/(1+(ξ+4*α^2))) * (D0*(1-((ξ+4*α^2)/2))^2))
αmean(1,1/3,0.0)

remake(prob, u0 = rand() * prob.u0)

=#

# Trajectories

α = 3^0.5/6
ξ = 1/3
tmax = 128
N = 2^14
dt = 1e-4
beta = 1

Random.seed!(123456)

dx = 2 * pi/N    
eta = 2 * dx
kappa =eta^beta
r = [-pi + i * dx for i in 0:N-1]
logker = SingularKernel(r, eta);
push!(r, pi)
σ_sq = 2 * sum(abs.(logker.Lk[2:end]).^2)

#setting environment
g2 = UnitaryWhiteNoise(div(length(logker.r), 2) + 1)
gmc = exp.(α .* realization(logker, g2) .- α^2 .* σ_sq)
gmcinterp = linear_interpolation(r, push!(gmc, gmc[1]), extrapolation_bc = Periodic())    
#hdf5()
#plt = plot(gmcinterp(r),r)

p=[gmcinterp]
prob_func = let p = p
    (prob, i, repeat) -> begin
        g2 = UnitaryWhiteNoise(div(length(logker.r), 2) + 1)
            gmc = exp.(α .* realization(logker, g2) .- α^2 .* σ_sq)
            gmcinterp = linear_interpolation(r, push!(gmc, gmc[1]), extrapolation_bc = Periodic())
        
            remake(prob, p = [gmcinterp])
        end
    end
    
function piecewisekernel(x, y, η, ξ)
    r=x-y
    @.((1. - normevans(r,η)^ξ) * (r<=1)) 
end

function f(du,u, p, t) 
    du[1]=0.
    du[2]=0.
end


function g(du, u, p, t)
    L_sq = [piecewisekernel(x, y, eta, ξ) for x in u, y in u]
    if isposdef(L_sq)
        L, U = cholesky(L_sq)
    else
        L = sqrt(L_sq)
    end
    gmcinterp = p[1]
    du[1,1] = L[1,1]*gmcinterp(u[1])
    du[1,2] = L[1,2]*gmcinterp(u[2])
    du[2,1] = L[2,1]*gmcinterp(u[1])
    du[2,2] = L[2,2]*gmcinterp(u[2])
    #kappa
    du[1,3] = (2*kappa)^0.5
    du[2,3] = (2*kappa)^0.5 
end

function condition(u, t, integrator) 
    abs(u[1]-u[2]) - 1
end

function condition2(u, t, integrator) 
    abs(u[1]-u[2]) < eta
end


affect1!(integrator) = terminate!(integrator)

function affect2!(integrator)  
    integrator.u[1]=0.5*(integrator.u[1]+integrator.u[2])+eta/2
    integrator.u[2]=0.5*(integrator.u[1]+integrator.u[2])-eta/2
end

cb1 = ContinuousCallback(condition, affect1!)
cb2 = SciMLBase.DiscreteCallback(condition2,affect2!)

cbs=CallbackSet(cb1,cb2)
u0 = [-2 ,-2 + eta]
tspan = (0.0, 50)
prob = SDEProblem(f, g, u0, tspan, p, noise_rate_prototype = zeros(2, 3))
#ensembleprob = EnsembleProblem(prob,output_func = (sol, i) -> (sol.t[end], false), prob_func = prob_func)
#global sol = solve(ensembleprob,EM(),EnsembleDistributed();dt=dt,callback=cb, trajectories = nindep)
sol = solve(prob, EM(),dt=dt,callback=cbs)
df = DataFrame(sol)
t = df.timestamp

x1 = df.value1
x2 = df.value2


d = Dict()
d["gmcinterp"]=gmcinterp
plot(gmcinterp(r),r)
plot!(t,[x1,x2])

sol2 =sol
sol3=sol
plot(gmcinterp(r),r)
hdf5()
plt=plot(sol)
Plots.hdf5plot_write(plt, "trajectories.hdf5")
plotlyjs() # Must first select some backend
pread = Plots.hdf5plot_read("trajectories.hdf5")



#plot
using Interpolations
using JLD2
file = jldopen("trajectories.jld2")

gmcinterp = file["gmcinterp"]
r = file["gmc_x"]
using ColorSchemes
using Plots;gaston(size=(800,650))
c=cgrad(:bluesreds,scale=:exp, rev = false, categorical = false);
default(titlefont = (20, "times"), legendfontsize = 18, guidefont = (12, :darkgreen), tickfont = (12, :orange), guide = "x", framestyle = :box, yminorgrid = true)

plt = plot(layout=(1,1));
a = 3/2
Plots.scalefontsizes(1/a)
plot(sin, thickness_scaling = 1)
plt=plot(-gmcinterp(r),r,
    seriestype = :path,
    line_z = gmcinterp(r),
    color = c,
    framestyle = :box,
    aspect_ratio = 1,
    widen = false,
    xlims = [-9,0],
    thickness_scaling = a,

    )
vline!([xlims(plt)[1], xlims(plt)[2]],lc = :black,lw = 5, label = "")
hline!([ylims(plt)[1], ylims(plt)[2]],lc = :black,lw = 5,label = "")


t=file["t56"]
x1=file["x5"]
x2=file["x6"]

plot!(t,x2,
     #seriestype = :path,
     #line_z = i-> get(ColorSchemes.rainbow,((df[!,"value1"][i]>0)+1)/2),   
     line_z = gmcinterp(x2),
     color = c,
     colorbar=:true,
)

y=randn(10)
