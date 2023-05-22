using DrWatson
@quickactivate "IntermittentFields";
include(srcdir("IntermittentFields.jl"))
using .IntermittentFields
using LinearAlgebra
using Statistics
using Plots
using DataFrames
using LsqFit
using StatsBase
using PyFormattedStrings
using LaTeXStrings
using Interpolations
using DifferentialEquations
using ColorSchemes

white_list_fields = ["N", "α", "p", "t", "ξ", "ncor", "nindep","np", "nproc", "λ","dt", "tmax"]
white_list_sde = ["N", "α", "p", "t", "ξ", "nindep","np", "nproc","dt", "tmax"]
white_list_sdesolver = ["sol","N", "α", "ξ", "nindep", "nproc","dt", "tmax"]

data_fields = collect_results(datadir("sims", "dispersion", "quenchedgmc"), white_list = white_list_fields)
data_sde = collect_results(datadir("sims", "dispersion", "frozenfield"), white_list = white_list_sde,rexclude=[r"diagGMC",r"noeta",r"expkernel",r"sqrt"])
data_sdesolver = collect_results(datadir("sims", "dispersion", "sdesolver"), white_list = white_list_sdesolver,rinclude=[r"solver2_"])

sort!(data_fields, :N)
sort!(data_sde, :N)   
sort!(data_sdesolver, :N)    


#transform!(data, :t => ByRow( t -> normalize(fit(Histogram, vec(t), 0:0.001:65), mode = :pdf)) => :tpdf)
transform!(data_fields, :t => ByRow(maximum) => :maxt)
transform!(data_fields, :t => ByRow(mean) => :avgt)
transform!(data_fields, :t => ByRow(t->mean(min.(t, 64))) => :avgt64)
transform!(data_fields, :t => ByRow(t->mean(min.(t, 32))) => :avgt32)
transform!(data_fields, :t => ByRow(t->mean(min.(t, 16))) => :avgt16)
transform!(data_fields, :t => ByRow(t->mean(min.(t, 8))) => :avgt8)

transform!(data_sde, :t => ByRow(maximum) => :maxt)
transform!(data_sde, :t => ByRow(mean) => :avgt)
transform!(data_sde, :t => ByRow(t->mean(min.(t, 128))) => :avgt128)
transform!(data_sde, :t => ByRow(t->mean(min.(t, 64))) => :avgt64)
transform!(data_sde, :t => ByRow(t->mean(min.(t, 32))) => :avgt32)
transform!(data_sde, :t => ByRow(t->mean(min.(t, 16))) => :avgt16)
transform!(data_sde, :t => ByRow(t->mean(min.(t, 8))) => :avgt8)

transform!(data_sdesolver, :sol =>  ByRow(sol->maximum(sol.u)) => :maxt)
transform!(data_sdesolver, :sol => ByRow(sol->var(min.(sol.u, 128))) => :var)
transform!(data_sdesolver, :sol => ByRow(sol->mean(min.(sol.u, 128))) => :avgt128)
transform!(data_sdesolver, :sol => ByRow(sol->mean(min.(sol.u, 64))) => :avgt64)
transform!(data_sdesolver, :sol => ByRow(sol->mean(min.(sol.u, 32))) => :avgt32)
transform!(data_sdesolver, :sol => ByRow(sol->mean(min.(sol.u, 16))) => :avgt16)
transform!(data_sdesolver, :sol => ByRow(sol->mean(min.(sol.u, 8))) => :avgt8)
#transform!(data_sdesolver, :sol => ByRow(sol -> normalize(fit(Histogram, vec(min.(sol.u, 128)), 0:0.1:128), mode = :pdf)) => :tpdf)

Plots.scalefontsizes(1/1.5)
Plots.default(#titlefont = (20, "times"),
            guidefont = (18,"serif-roman"),
            tickfont = (18, :black),
            framestyle = :box,
            fontfamily="serif-roman",
            #yminorgrid = true,   
)

#τ X α FIG 2c ----------------------------------------------------------
begin
    using Plots.PlotMeasures
    gr(size=(800,600))
    mylabels = [ξ -> f"Fields", ξ -> f"SDE"]
    mylabels_theoretical = [L"\xi = 1/3", L"\xi = 2/3"] 
    mymarker = [:circle, :xcross,];
    inside_color=[:white,:white]

    
    mystrokewidth = [2, 2]
    mycolors = [:firebrick, :royalblue]
    ξs = [1/3, 2/3];
    myαs = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6] 
    αcrit = [sqrt((1 - ξ)/4) for ξ in ξs]
    mydata = [data_fields, data_sdesolver];
    theoretical(ξ,α) = 1/((1-ξ-4*α^2)*(2-ξ-4*α^2))

    p = plot(layout=(1,1));
    #prediction
    for (i, ξ) in enumerate(ξs) 
        plot!(0:0.0001:αcrit[i], α->theoretical(ξ,α),
        #ls = :dash,
        lc = mycolors[i],
        lw = 2,
        label = mylabels_theoretical[i],
        )
        vline!([αcrit[i]], ls = :dash, lc = mycolors[i], lw = 2, label = "")
    end 
    #@show p
    for (k, data) in enumerate(mydata), (j, ξ) in enumerate(ξs)
        df = data[(data.N .== 2^14) .& (data.ξ .== ξ) .& (data.dt .== 1e-4), :];
        sp = 1
        αs = []
        τ = []
        for i in 1:nrow(df)
            if df.α[i] in myαs
                αs = append!(αs, df.α[i])
                τ = push!(τ, df.avgt64[i])
            end
        end
        # 

        if (data==data_fields)
            global inside_color[j]= :white
        else
            inside_color[j]=mycolors[j]
        end
        #filter!((sqrt(3)/0.6), αs)

        #println(αs)
        sort!(αs)
        sort!(τ)
        plot!(αs,τ, 
            sp = sp,
            lw =1/4,
            lc=mycolors[j],
            #seriestype = :scatter,
            label =  mylabels[k](ξ),#"Fields-α = $α",
            markershape = mymarker[k],
            markersize = 6,
            markeralpha = 1,
            markercolor = inside_color[j], 
            markerstrokewidth = mystrokewidth[k],
            #markerstrokealpha = 0.1,
            markerstrokecolor = mycolors[j],
            #markerstrokestyle = :dash,
            xlabel = L"\gamma",
            xlims = (0,0.6),
            ylabel =  L"\tau_1",
            ylims = (0,70),
            yticks =[10, 20, 30, 40, 50, 60, 70] ,
            widen = true,
            thickness_scaling = 1.5,
            legend = :none,
            bottom_margin = 0mm,
            left_margin = 0mm,
            annotations = (0.15, 55, text("Mean Field", 12, "serif-roman")),
            #aspect_ratio=1/100
        )
        plot!([0.23,0.27],[55,53], arrow = true,color=:black,linewidth=1,label="")
        plot!([0.23,0.39],[55,65], arrow = true,color=:black,linewidth=1,label="")

    end
    
    @show p
end

#<τ> x η fig2a, fig2b---------------------------------------------------
begin
    using Plots.PlotMeasures
    gr(size = (800, 600))
    mylabels = ["Fields" "Fields"; "SDE" "SDE"]
    mylabels_x = ["" ""; L"\eta" L"\eta"]
    mylabels_y = [L"\left<\tau\right>" ""; L"\left<\tau\right>" ""]
    mylabels_tmax = [L"\tau\wedge 8", L"\tau\wedge 16",L"\tau\wedge 32",L"\tau\wedge 64"]
    mylabels_theoretical = [L"\xi = 1/3"  L"\xi = 2/3"; L"\xi = 1/3"  L"\xi = 2/3"] 
    mymarker = [:utriangle, :square, :diamond, :circle ];
    mystrokewidth = [1, 2]
    mycolors = [:orange, :teal]
    mysubplots = [1 2 ; 3 4]
    ξ = 2/3;
    αs = [0.0, 0.6];
    mydata = [data_fields, data_sdesolver];
    D0 = 1
    αmean(λ, ξ, α) = (λ^2)/(2 * ((2/(2-(ξ+4*α^2)))*(1-(ξ+4*α^2)^2)/(1+(ξ+4*α^2))) * (D0*(1-((ξ+4*α^2)/2))^2))
    mytheoretical_hline = αmean(1,ξ,0)

    #prediction
    for (k, data) in enumerate(mydata)
        p = plot(layout=(1,1));

        for (l, α) in enumerate(αs)
            df = data[(data.α .== α) .& (data.ξ .== ξ) .& (data.dt .== 1e-4) , :];
            
            #sp = mysubplots[k,j]
            hline!([mytheoretical_hline], ls = :dash, lc = :black, label = "")
            #hline!(sp = sp, [64], ls = :dashdot, lc = :black, label = "")
            ηs = []
            t64 = []
            t32 = []
            t16 = []
            t8 = []
            for i in 1:nrow(df)
                ηs = append!(ηs, 8pi/df.N[i])
                t64 = push!(t64, df.avgt64[i])
                t32 = push!(t32, df.avgt32[i])
                t16 = push!(t16, df.avgt16[i])
                t8 = push!(t8, df.avgt8[i])
            end
            for (m,t) in enumerate([t8, t16, t32, t64])
                plot!(ηs,t, 
                    sp = 1,
                    seriestype = :scatter,
                    label = mylabels_tmax[m], #L"\alpha = " * f"{α:0.1f}", #mylabels[k, j],#"Fields-α = $α",
                    markershape = mymarker[m],
                    markersize = 7,
                    #markeralpha = 1,
                    markercolor = mycolors[l], #get(ColorSchemes.rainbow, j/length(αs)),
                    #markerstrokewidth =  mystrokewidth[k],
                    ##markerstrokealpha = 0.1,
                    markerstrokecolor = :black,#mycolors[l],
                    #markerstrokestyle = :dash,
                    xlabel =  L"\eta",
                    xscale = :log10 ,
                    yscale = :log10,
                    xlims = (1e-4,1),
                    ylims = (1,1e2),
                    ylabel =  L"\tau_1" ,#mylabels_y[k, j],
                    #yticks =[10, 20, 30, 40, 50] ,
                    #tickfontsize = 5,
                    #xguidefontsize = 10,
                    #yguidefontsize = 10,
                    #widen = true,
                    thickness_scaling = 1.5,
                    #legend = :topleft,
                    annotations=[(5e-4, 2e1, text("Increasing  "*L"T_{max}", 12, "serif-roman",rotation = 90)),],
                    top_margin = 4mm,
                    right_margin = 4mm,
                    legend=:none
                
                )
            end
            
        end
        #arrow
        x, y = 1e-3, 5    # vector origin
        u, v = 0, 50      # vector
        xl, yl = (1e-4, 1), (1, 100)  # loglog plot limits (have to be modified to max range)
        # plot arrowheads as rotated unicode characters:
        lminmax = (minimum([xl[1];yl[1]]), maximum([xl[2];yl[2]]))
        plot!([x,x+u], [y,y+v], lc=:black,scale=:log10)
        #lxl, lyl = log10.(Plots.xlims(p)),  log10.(Plots.ylims(p))
        #sz =  pp.attr[:size]
        #ar = (lxl[2]-lxl[1])/(lyl[2]-lyl[1])*sz[2]/sz[1]
        #r = atan(ar*log10.(1+v/y), log10.(1+u/x)) * 180/pi    # rotation angle
        annotate!(x+u,y+v, text('\u25B2', 9, :black, rotation=0))

        if k==1
            savefig(p,"fig2a.pdf")
        else
            savefig(p,"fig2b(callback).pdf")
        end
    end
    #@show p
    #gr(size=(200,150))
    ##p1 = plot(map(plot, [p[1], p[3]])..., layout=(2,1))
    #p2 = plot(map(plot, [p[2], p[4]])..., layout=(2,1))
    #g1 = plot(p2[1])
    #g2 = plot(p2[2])
    #savefig(g2,"fig2b.pdf")
    #@show g1
end
@show p

# heatmap fig2d---------------------------------------------------------
begin
    
    using Dierckx
    using Plots.PlotMeasures
    
    gr(size=(800,600))
    df = data_sdesolver[(data_sdesolver.N .==2^14) .& (data_sdesolver.dt .== 1e-4),:]
    
    ξ0s = [2/9, 3/9, 4/9, 5/9, 6/9, 7/9, 8/9, 9/9]
    α0s = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
    f(ξ,α)=last(df[(df.α .== α) .& (df.ξ .== ξ), :].avgt64)/64
    #f(ξ,α)=sqrt(last(df[(df.α .== α) .& (df.ξ .== ξ), :].var))/128
    #g(ξ,α)=last(df[(df.α .== α) .& (df.ξ .== ξ), :].maxt)
    
    ztmp = f.(ξ0s,α0s')
    itp = Spline2D(ξ0s, α0s, ztmp;  kx = 3, ky = 3, s = 0.0);
    x = collect(range(extrema(ξ0s)..., length = 1000))
    y = collect(range(extrema(α0s )..., length = 1000))
    z = [itp(x,y) for x in x, y in y]
    
    c = cgrad(:jet,scale=:exp, categorical = false,rev = true);
    clim = (0.0,1);

    # HEATMAP
    p = heatmap(x, y, z',xlims=(2/9,1), color = c, clim = clim,xlabel = L"\xi",ylabel=L"\gamma", thickness_scaling = 1.5,
    left_margin = -9mm,
    bottom_margin = -5mm,
    right_margin = 6mm,
    colorbar = true,
    annotations = (0.5, 0.25, text("Mean Field", 12, "serif-roman")),
    )

    # GRID POINTS
    scatter!(p,[x for _ in α0s for x in ξ0s], [y for y in α0s for _ in ξ0s];
    markercolor = :white,
    markershape = :circle,
    ms = 2,
    label = "",
    clim = clim,
    c = c,
    markerstrokecolor =:white,
    widen = false,
    thickness_scaling = 1.5,
    )

    # MEAN FIELD
    αc(ξ) = sqrt((1-ξ)/4)
    ξs = 2/9:1e-4:3/3
    plot!(ξs, αc.(ξs),
    lw = 4,
    linestyle = :dash, 
    lc = :black ,
    label="",
    widen=false,
    ylims=[0,0.60001],
    
    )  

    #ARROW
    x, y = 0.5, 0.27    # vector origin
    u, v = 0., 0.06     # vector
    xl, yl = (0, 0), (1, 0.5)  #plot limits (have to be modified to max range)
    lminmax = (minimum([xl[1]; yl[1]]), maximum([xl[2]; yl[2]]))
    plot!([x,x+u], [y,y+v], lc = :black, label = "")
    annotate!(x+u,y+v, text('\u25B2', 9, :black, rotation=0))
end 
@show p

#FIG 3c, 3d ------------------------------------------------------------
begin 
    using Plots.PlotMeasures
    gr(size = (800, 600))
    
    #dt = 1/10
    #transform!(data_sdesolver, :sol => ByRow(sol -> normalize(fit(Histogram, vec(min.(sol.u, 128)), 0:dt:128), mode = :pdf)) => :tpdf)
    
    dt_log = 1/2
    dt_lin = 1/10
    transform!(data_sdesolver, :sol => ByRow(sol -> normalize(fit(Histogram, vec(min.(sol.u, 128)), 0:dt_log:128), mode = :pdf)) => :logpdf)
    transform!(data_sdesolver, :sol => ByRow(sol -> normalize(fit(Histogram, vec(min.(sol.u, 128)), 0:dt_lin:128), mode = :pdf)) => :linpdf)
    
    df = data_sdesolver[(data_sdesolver.N .== 2^14) .& (data_sdesolver.α .== 3^0.5/6) .& (data_sdesolver.ξ .== 1/3) .& (data_sdesolver.dt .== 1e-4),:]
    df = data_sdesolver[(data_sdesolver.N .== 2^14) .& (data_sdesolver.α .== 0) .& (data_sdesolver.ξ .== 2/3).& (data_sdesolver.dt .== 1e-4),:]
    
    x_lin = Array(range(0, stop = 128, length = Int64(128/dt_lin)))
    x_log = Array(range(0, stop = 128, length = Int64(128/dt_log)))
   
    plot([x_lin, x_log],  [df.linpdf[1].weights, df.logpdf[1].weights],
        seriestype = [:path, :path],
        inset_subplots = bbox(0.55, 0.058, 0.4, 0.3),
        labels = ["" ""],
        xlims = [(0, 10)  (0, 64)],
        ylims = [(0, 0.4) (1e-4, 1e0)], 
        #xticks=[ 0:0.5:10  0:0.1:10],
        #xscale = [:linear, :log],
        yscale = [:identity  :log10],
        lw =[4  4],
        lc = [:royalblue  :royalblue],
        xlabel = [L"T_{1}"  ""],
        ylabel = ["PDF"  ""],
        thickness_scaling = 1.5,
        bottom_margin = 0mm,
        left_margin = 0mm,
        right_margin = 3mm,
        #annotations=[(1, 1, text("Increasing  "*L"T_{max}", 12, "serif-roman",rotation = 90)),],
    
    );
    #plot!(sp = 2, yticks = nothing)
    annotate!(sp = 1, 7, 0.1, text(L"\propto \tau_{1}^{-1}e^{-T_1/\tau_1}", 12, "serif-roman", rotation = 0));
    plot!(sp = 2, xtickfontsize = 9, ytickfontsize = 9);
    
    #Exponential
    #y_lin = 1/(df.avgt128[1])*exp.(-x_lin/df.avgt128[1]);
    #y_log = 1/(df.avgt128[1])*exp.(-x_log/df.avgt128[1]);
    dt = 1/100
    x = Array(range(0, stop = 128, length = Int64(ceil(128/dt))))
    y = 1/(df.avgt128[1])*exp.(-x/df.avgt128[1]);
    
    plot!(x, [y, y],
        sp=[1 2],
        lw = [4 4], 
        linestyle = [:dash :dash],
        lc = [:black :black], 
        label=["" ""],    
    )

end

# FIG 1a ---------------------------------------------------------------
begin
    gr(size = (800, 600))
    
    #loading r0s
    df = collect_results(datadir("sims", "field_statistics"), rinclude=[r"structurefunc"])
    r0s = df.l[1]
    
    #lognormal curve
    ξ = 2/3
    α =0.2
    ps = 0:0.01:6
    ζs = @. (ps*(ξ/2 + α^2) -α^2*ps^2/2)
    
    #monofractal
    mono = ps .* ξ/2
    
    plot(ps, [ζs, mono],
        sp = 1,
        ls = [:solid :dash],
        lc = [:black :black],
        lw = [2 1],
        lab = ["" ""],
        xlabel = L"p",
        ylabel = L"\zeta(p)",
        inset_subplots = bbox(0.30, 0.058, 0.2*1.5, 0.15*1.5),
        thickness_scaling = 1.5
    );
    
    # subplot correlation function
    r = collect(-pi:1e-4:pi)
    cr = @. (1 - abs(r)^ξ)*(abs(r)<=1)
    cr0s = @. (1 - abs(r0s)^ξ)*(abs(r0s)<=1)
    #eta = 0
    #ker = CovarianceKernel(r, eta, CovarianceCorrelation(piecewisekernel), ξ, false);
    
    #markers
    x = collect(1:1:6);
    m = @. (x*(ξ/2 + α^2) -α^2*x^2/2);
    mc = get(cgrad(:rainbow, length(x), rev = false, categorical = true), 2)
    
    #plot(r,[ker.Lr_sq,cr])
    plot!(sp = 2, [r, r0s], [cr, cr0s],
        st = [:path :scatter],
        xlims = [0, 1],
        ylims = [0, 1],
        #xscale = :log10,
        #yscale = :log10,
        lw = 1,
        #ls = [:solid],
        #linealpha = 1,
        lc = :black,
        mc = [get(cgrad(:rainbow, length(x), rev = false, categorical = true), 2/length(x)) for _ in eachindex(x)],
        lab = "",
        ylabel = L"C_\xi(r)",
        xlabel = L"r",
        xtickfontsize = 9, 
        ytickfontsize = 9,
        xticks = [0, 1],
        yticks = [0, 1],
        widen = false,
        xguidefontsize = 10,
        yguidefontsize = 10
    );
    
    
    c = cgrad(:rainbow, length(x), rev = false, categorical = true);
    plot!(x, m,
        seriestype = :scatter,
        cbar = false,
        c = [get(cgrad(:rainbow, length(x), rev = false, categorical = true), i/length(x)) for i in eachindex(x)],
        ms = 6,
        label = "",
    )
    
    annotate!(sp = 1, 4.5, 0.8, text(L"p \,\, \left(\frac{\xi}{2} + \gamma^2 \right) - \frac{\gamma^2 p^2}{2}", 12, "serif-roman", rotation = 0));
    annotate!(sp = 1, 5, 1.9, text(L"p \frac{\xi}{2}", 12, "serif-roman", rotation = 0));
    
    #arrows
    plot!(sp = 1, [4.5, 4.5],[0.95, 1.2], arrow = true, label = "", c = :black, lw =1)
    plot!(sp = 1, [5.2, 5.5],[1.9, 1.9], arrow = true, label = "", c = :black, lw =1)
    
end

#FIG 1b
begin
    using Plots.Measures
    gr(size = (800, 600))

    df=collect_results(datadir("sims", "field_statistics"), rinclude=[r"structurefunc"])
    sp = collect(df.sp_abs[1])
    r0s = df.l[1]
    ξ = df.ξ[1]
    α = df.α[1]
    eta = df.eta[1]
    ps = df.ps[1]
    samples = df.samples[1]
    ydata=[r0s.^(-ξ).*sp[i].^(2/ps[i]) for i in eachindex(ps)]
    
    #lognormal prediction
    ζs = @. (ps*(ξ/2 + α^2) -α^2*ps^2/2)
    theo = ((2ζs./ps).-ξ)
    v=[linear_interpolation(r0s,ydata[i]) for i in eachindex(ydata)]
    
    p = plot(r0s, ydata,#[v[i](1)*r0s.^th[i] for i in eachindex(v)],
            marker = :circle,        
            axis = :log,
            ms = 6,
            lw =2,
            #c = :black,
            #linealpha=0.3, 
            xlims = [1e-6, 10],
            ylims = [1e-1, 1e2],
            xticks = [1e-6, 1e-4, 1e-2, 1],
            markercolor= [get(cgrad(:rainbow, length(ps),rev=false,categorical=true), i/length(ps)) for i in eachindex(ps)]',
            lc= [get(cgrad(:rainbow, length(ps),rev=false,categorical=true), i/length(ps)) for i in eachindex(ps)]',
            thickness_scaling = 1.5,
            label = "",
            xlabel = L"r",
            ylabel = L"r^{-\xi}\left(\mathrm{E}|\delta u|^p \right)^{2/p}", 
            left_margin = -10mm,
            top_margin = 4mm,

    );
    
    annotate!(0.01, 0.5, text("Inertial range", 12,"serif-roman"));
    annotate!(3e-6, 7e1, text(L"r<\eta", 12,"serif-roman"));
    annotate!(3e0, 7e1, text(L"1<r", 12,"serif-roman")); 
    a=[1.2, 2, 2.8, 3.7, 4.7, 6]
    plot!(r0s[3:end-2],[a[i]*r0s[3:end-2].^theo[i] for i in eachindex(v)], lab = "", lc = :black, lw=3, la=0.4)
    vspan!([xlims(p)[1], eta ], color = :grey, alpha = 0.2, labels = "");
    vspan!([1,xlims(p)[2] ], color = :grey, alpha = 0.2, labels = "")

    
    #fit
    #@. model(r,c)= c[1]*r^c[2]
    #fit = curve_fit(model,r0s[6:end-2], ydata[6][6:end-2],[1.,1.])
    #plot!(r0s, @. 1.2*fit.param[1]*r0s^fit.param[2])
    

end

#Fig 1c -----------------------------------------------------------------
begin
    gr(size = (800, 600))
    
    df = collect_results(datadir("sims", "field_statistics"), rinclude = [r"realizations"])
    r = df[!,"r"][1]
    u = df[!,"fields"][1]
    αs = df[!,"αs"][1]
    
    x0 = 4.5
    dx = 8.5
    p = plot(xlims = [-pi, pi], ylims = [-4, 40], thickness_scaling = 1.5);
    for (i, α) in enumerate(αs)
        plot!(r, u[i] .+(i-1)*8.5,
            lab = "",
            lw = 2,
            c = get(cgrad(:ice, 2*length(αs),rev = true,categorical=true), (2i)/(2*length(αs))),
            xlabel = L"x",
            ylabel = L"u_{\eta}^{\xi,\gamma}",
        )
        #annotate!(-2, x0+ (i-1)*dx, text(L"\gamma="*f"{α:0.1f}", 12, "serif-roman"));
    end
    annotate!(-2, x0, text(L"\gamma=0", 12, "serif-roman"));
    annotate!(-2, x0 + 1dx, text(L"\gamma=0.2", 12, "serif-roman"));
    annotate!(-2, x0 + 2dx, text(L"\gamma=0.4", 12, "serif-roman"));
    annotate!(-2, x0 + 3dx, text(L"\gamma=0.6", 12, "serif-roman"));
    @show p
end

#FIG 1d -----------------------------------------------------------------
begin
    gr( size = (800, 600))
    using Plots.Measures
    
    ζ(p, ξ, α) = p*(ξ/2 + α^2) - α^2*p^2/2
    ξk(p, ξ, α) = 2*(ζ(p, ξ, α)-1)/p
    
    ps = 0.1:1e-2:1000
    maxξ(ξ, α) = maximum(ξk.(ps, ξ, α))
    #mm(3/3,0)
    αs = 0:1e-3:0.6
    plot(αs, [maxξ.(1, αs), maxξ.(2/3, αs), maxξ.(1/3, αs)],
        inset_subplots = bbox(0.54, 0.05, 0.4, 0.3),
        sp = 1 ,
        xlims = [0, 0.6], 
        ylims = [0, 1], 
        yticks = ([0, 1/3, 2/3, 1], ["0", L"\frac{1}{3}",  L"\frac{2}{3}", "1"]),
        xticks= ([0, 0.2, 0.4, 0.6],["0", "0.2", "0.4", "0.6"]),
        widen = false,
        lab = "",
        lc = :royalblue,
        lw = 2,
        thickness_scaling = 1.5,
        xlabel = L"\gamma",
        ylabel = L"\xi_K",
        right_margin = 3mm,
        left_margin = -3mm
    );
    
    ps = 0.1:1e-2:1000
    ξ = 2/3
    α = 0.1
    plot!(sp = 2, ps, ξk.(ps, ξ, α),
        xlims = [0, 70],
        ylims = [0, 1.05ξ],
        widen = false,
        lab = "",
        lc = :black,
        lw = 2,
        xticks = ([0, 2^0.5/α, 70],["0", L"p_K=\sqrt{2}/\gamma", L"p"]),
        yticks = ([maximum(ξk.(ps, ξ, α)), ξ],[L"\xi_K", L"\xi"]),
        xtickfontsize = 9, 
        ytickfontsize = 9,
    
    );
    hline!(sp = 2, [ξ], lw = 2, ls = :dash, lc = :black, lab = "");
    plot!(sp = 2, [2^0.5/α, 2^0.5/α], [0, maximum(ξk.(ps, ξ, α))], lc = :royalblue, lw = 2, ls = :dash, lab = "")
    plot!(sp = 2, [0, 2^0.5/α], [maximum(ξk.(ps, ξ, α)), maximum(ξk.(ps, ξ, α))], lc = :royalblue, lw = 2, linealpha = 1,ls = :dash, lab = "")
    scatter!(sp = 2, [2^0.5/α], [maximum(ξk.(ps, ξ, α))], c = :royalblue,lab = "")
    annotate!(sp = 2, 60, 0.3, text(L"2\frac{\zeta(p)-1}{p}", 9, "serif-roman"))
    plot!(sp = 2,[48, 40],[0.3, 0.3] ,arrow = true, lc = :black, lab = "")
end


###################################
    

