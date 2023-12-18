function plot_multiple_AEW(onlines; save=false, filename="figs/LDH_AEW_raw.png")
    ps = [scatter(online[:,1]./60, online[:,2]) for online in onlines]
    ps_tot = plot([p for p in ps]...,size=(3000,3000))
    display(ps_tot)
    if save
        savefig(ps_tot, filename)
    end
end

function plot_intensity(online; save=false, filename="figs/intensity.png")
    p = plot(online[:,1]./60, online[:,3], size=(500,350),
        label = "",
        legendfontsize = 7,
        titlelocation = :left,
        bottom_margin=10Plots.px,
        left_margin=10Plots.px,
        tickfontsize = 10,
        xlabelfontsize = 10,
        ylabelfontsize = 10,
        grid = false,
        framestyle = :box,
        #legend=[:topleft false false :topleft false false],
        ylabel = "Intensity",
        xlabel = "Time in hours",
        fillrange = online[1,3],
        fillalpha = 0.35,
        c = 4,
        linewidth =0,
        ylim=(2E5,3E5),
        xlim=(-0.03,online[end,1]./60),
        title = "Intensity of LDH 7")
    plot!(online[:,1]./60, online[:,3], label="measured intensity", c=1, linewidth=2.5)
    scatter!([online[1,1]./60], [online[1,3]], label="Initial F at 100% soluble P", c=1, markersize=7)
    annotate!(1.75, 2.7E5, text("Insoluble \nAggregates", 4, :right, 15))
    annotate!(0.65, 2.15E5, text("Soluble \nprotein fraction", :black, :center, 15))
    #annotate!(1.75, 2.7E5, text("Aggregates", 4, :right, 15))
    #plot!(online[[1,end],1]./60, online[[1,1],3], label = "estimated intensity")
    display(p)
    if save
        savefig(p, filename)
    end
end

function plot_intensity_2(online; save=false, filename="figs/intensity.png")
    p = plot(online[:,1]./60, online[:,3], size=(500,350),
        label = "",
        legendfontsize = 7,
        titlelocation = :left,
        bottom_margin=10Plots.px,
        left_margin=10Plots.px,
        tickfontsize = 10,
        xlabelfontsize = 10,
        ylabelfontsize = 10,
        grid = false,
        framestyle = :box,
        #legend=[:topleft false false :topleft false false],
        ylabel = "Intensity",
        xlabel = "Time in hours",
        fillrange = online[1,3],
        fillalpha = 0.35,
        c = 4,
        linewidth =0,
        #ylim=(2E5,3E5),
        xlim=(-0.03,online[end,1]./60),
        title = "(A) Intensity of GalOx 3")
    plot!(online[:,1]./60, online[:,3], label="measured intensity", c=1, linewidth=2.5)
    scatter!([online[1,1]./60], [online[1,3]], label="Initial F at 100% soluble P", c=1, markersize=7)
    annotate!(1.75, 1.9E5, text("Cofactor \nquenching", 4, :center, 15))
    #annotate!(0.65, 2.15E5, text("Soluble \nprotein fraction", :black, :center, 15))
    #annotate!(1.75, 2.7E5, text("Aggregates", 4, :right, 15))
    #plot!(online[[1,end],1]./60, online[[1,1],3], label = "estimated intensity")
    display(p)
    if save
        savefig(p, filename)
    end
    return p
end

function plot_multiple_dAEW_adapted(onlines; save=false, filename="figs/LDH_dAEW_adapted.png")
    ps = []
    for online in onlines
        t = online[:,1]./60
        push!(ps, scatter(t,online.dAEW_adapted))
    end
    ps_tot = plot([p for p in ps]...,size=(3000,3000))
    display(ps_tot)
    if save
        savefig(ps_tot, filename)
    end
end

function plot_multiple_diff_aew(onlines; save=false, filename="figs/LDH_diff_AEW.png")
    ps = []
    for online in onlines
        t = online[:,1]./60
        p=scatter(t,online.diff_AEW_raw, markersize = 6, markeralpha=0.4, markerstrokewidth = 0.2, label="finite-difference")
        scatter!(t,online.diff_AEW_sgol, markersize = 6, markerstrokewidth = 0.2, label="Savitzky-Golay")
        scatter!(t,online.diff_AEW_loess, markersize = 6, markerstrokewidth = 0.2, label="Loess")
        push!(ps, p)
    end
    ps_tot = plot([p for p in ps]...,size=(3000,3000))
    display(ps_tot)
    if save
        savefig(ps_tot, filename)
    end
end

function plot_multiple_integrals_AEW(onlines; save=false, filename="figs/LDH_compare_integrals.png")
    ps = []; RMSE_raw = []; RMSE_sgol = []; RMSE_loess = []
    for online in onlines
        p= scatter(online[2:end,1]./60, -cumsum(online.diff_AEW_raw[2:end].*diff(online[:,1]./60)), markersize = 4, markerstrokewidth = 0.2, label="finite-difference")
        scatter!(online[2:end,1]./60, -cumsum(online.diff_AEW_sgol[2:end].*diff(online[:,1]./60)), markersize = 4, markerstrokewidth = 0.2, label="Savitzky-Golay")
        scatter!(online[2:end,1]./60, -cumsum(online.diff_AEW_loess[2:end].*diff(online[:,1]./60)), markersize = 4, markerstrokewidth = 0.2, label="Loess")
        plot!(online[:,1]./60, online.dAEW_adapted, label="measured", color="red", linewidth=3)
        push!(ps, p)
        push!(RMSE_raw, (sqrt(mean((online.dAEW_adapted[2:end] .- (-cumsum(online.diff_AEW_raw[2:end].*diff(online[:,1]./60)))).^2))) ./ ((maximum(online.diff_AEW_raw[2:end]) - minimum(online.diff_AEW_raw[2:end]))) )
        push!(RMSE_sgol, (sqrt(mean((online.dAEW_adapted[2:end] .- (-cumsum(online.diff_AEW_sgol[2:end].*diff(online[:,1]./60)))).^2))) ./ ((maximum(online.diff_AEW_sgol[2:end]) - minimum(online.diff_AEW_sgol[2:end]))) )
        push!(RMSE_loess, (sqrt(mean((online.dAEW_adapted[2:end] .- (-cumsum(online.diff_AEW_loess[2:end].*diff(online[:,1]./60)))).^2))) ./ ((maximum(online.diff_AEW_loess[2:end]) - minimum(online.diff_AEW_loess[2:end]))) )
    end
    ps_tot = plot([p for p in ps]...,size=(3000,3000))
    display(ps_tot)
    println("sum RMSE raw: ", mean(RMSE_raw))
    println("sum RMSE sgol: ", mean(RMSE_sgol))
    println("sum RMSE loess: ", mean(RMSE_loess))
    # h = histogram(RMSE_raw, label="finite-difference", bins=0:0.05:0.5, alpha=0.5)
    # h = histogram!(RMSE_sgol, label="Savitzky-Golay", bins=0:0.05:0.5, alpha=0.5)
    # h = histogram!(RMSE_loess, label="Loess", bins=0:0.05:0.5, alpha=0.5)
    #display(h)
    if save
        savefig(ps_tot, filename)
    end
end

function plot_AEW_vs_dAEW(onlines; save=false, filename="figs/plot_AEW_vs_dAEW.png")
    ps = []; ps2=[];
    for (i, online) in enumerate(onlines)
        t = online[:,1]./60
        push!(ps2, scatter(t, online[:,2], markersize = 3, markerstrokewidth = 0.2, label="", title="LDH $(i)"))
        if i == 1
            ylabel!("AEW in nm")
        end
        p=scatter(t,online.diff_AEW_raw, markersize = 3, markeralpha=0.4, markerstrokewidth = 0.2, label="Finite Difference",
            xlabel="Time in hours")
        if i == 1
            ylabel!(L"$-\frac{d}{dt}$(AEW) in nm/hour")
        end
        scatter!(t,online.diff_AEW_sgol, markersize = 3, markerstrokewidth = 0.2, markershape = :diamond, label="Savitzky-Golay")
        scatter!(t,online.diff_AEW_loess, markersize = 3, markerstrokewidth = 0.2, markershape = :xcross, label="Loess")
        push!(ps, p)
    end
    ps_tot = plot([p for p in vcat(ps2,ps)]..., size=(1000,550), layout=(2,length(onlines)),
        #legendfontsize = 10,
        titlelocation = :left,
        bottom_margin=10Plots.px,
        left_margin=10Plots.px,
        tickfontsize = 10,
        xlabelfontsize = 10,
        ylabelfontsize = 10,
        grid = false,
        framestyle = :box)
    display(ps_tot)
    if save
        savefig(ps_tot, filename)
    end
    return ps_tot
end

function plot_AEW_vs_dAEW_b(onlines, t_cop; save=false, filename="figs/plot_AEW_vs_dAEW.png")
    ps = []; ps2=[];
    for (i, online) in enumerate(onlines)
        t = online[:,1]./60
        p = plot(t, online[:,2], label="", title="Process $(i)",ylabel = "AEW in nm")
        #vline!([t_cop[i]/60], label="t_cop", color="red")
        #push!(ps2, p2)

        subplot = twinx()
        p = scatter!(subplot,t,online.diff_AEW_raw, markersize = 3, markeralpha=0.4, 
            markerstrokewidth = 0.2, label="Finite Difference",
            xlabel="Time in hours", ylim=(0,10), ylabel = "dAEW in nm/hour")
        scatter!(subplot,t,online.diff_AEW_sgol, markersize = 3, markerstrokewidth = 0.2, markershape = :diamond, label="Savitzky-Golay")
        scatter!(subplot,t,online.diff_AEW_loess, markersize = 3, markerstrokewidth = 0.2, markershape = :xcross, label="Loess")
        vline!(subplot,[t_cop[i]/60], label="t_cop", color="red")
        push!(ps, p)
    end
    ps_tot = plot([p for p in ps]..., size=(1000,550), layout=(2,Int64(length(onlines)/2)),
        legendfontsize = 5,
        titlelocation = :left,
        bottom_margin=10Plots.px,
        #left_margin=10Plots.px,
        tickfontsize = 10,
        xlabelfontsize = 10,
        ylabelfontsize = 10,
        grid = false,
        framestyle = :box,
        )
    # ps_tot2 = plot([p for p in ps2]..., size=(1000,550), layout=(2,Int64(length(onlines)/2)),
    #     legendfontsize = 7,
    #     titlelocation = :left,
    #     bottom_margin=10Plots.px,
    #     left_margin=10Plots.px,
    #     tickfontsize = 10,
    #     xlabelfontsize = 10,
    #     ylabelfontsize = 10,
    #     grid = false,
    #     framestyle = :box,
    #     #legend=[false true false false false false]
    #     )
    display(ps_tot)
    # display(ps_tot2)
    if save
        savefig(ps_tot, filename)
    end
end


function plot_AEW_vs_dAEW_c(onlines, t_cop; save=false, filename="figs/Galox_dAEW.pdf", filename2="figs/Galox_AEW.pdf")
    ps = []; ps2=[];
    for (i, online) in enumerate(onlines)
        t = online[:,1]./60
        p2 = scatter(t, online[:,2], label="", title="GalOx $(i)",ylabel = "AEW in nm", markersize = 3,
        markerstrokewidth = 0.7)
        vline!([t_cop[i]/60], label="copper addition", color="red")
        push!(ps2, p2)

        p = scatter(t,online.diff_AEW_raw, markersize = 3, markeralpha=0.4, title="GalOx $(i)",
            markerstrokewidth = 0.2, label="Finite Difference",
            xlabel="Time in hours", ylim=(0,10), ylabel = "dAEW in nm/hour")
        scatter!(t,online.diff_AEW_sgol, markersize = 3, markerstrokewidth = 0.2, markershape = :diamond, label="Savitzky-Golay")
        scatter!(t,online.diff_AEW_loess, markersize = 3, markerstrokewidth = 0.2, markershape = :xcross, label="Loess")
        vline!([t_cop[i]/60], label="copper addition", color="red")
        push!(ps, p)
    end
    ps_tot = plot([p for p in ps]..., size=(1000,550), layout=(2,Int64(length(onlines)/2)),
        legendfontsize = 7,
        titlelocation = :left,
        bottom_margin=10Plots.px,
        left_margin=10Plots.px,
        tickfontsize = 10,
        xlabelfontsize = 10,
        ylabelfontsize = 10,
        grid = false,
        framestyle = :box,
        legend=[:topleft false false :topleft false false],
        ylabel = ["dAEW in nm/hour" "" "" "dAEW in nm/hour" "" ""],
        xlabel = ["" "" "" "Time in hours" "Time in hours" "Time in hours"],
        )
    ps_tot2 = plot([p for p in ps2]..., size=(1000,550), layout=(2,Int64(length(onlines)/2)),
        legendfontsize = 7,
        titlelocation = :left,
        bottom_margin=10Plots.px,
        left_margin=10Plots.px,
        tickfontsize = 10,
        xlabelfontsize = 10,
        ylabelfontsize = 10,
        grid = false,
        framestyle = :box,
        legend=[:topleft false false :topleft false false],
        ylabel = ["AEW in nm" "" "" "AEW in nm" "" ""],
        xlabel = ["" "" "" "Time in hours" "Time in hours" "Time in hours"],
        )
    display(ps_tot)
    display(ps_tot2)
    if save
        savefig(ps_tot, filename)
        savefig(ps_tot2, filename2)
    end
end

function plot_AEW_vs_dAEW_d(onlines, t_cop; save=false, filename="figs/Galox_dAEW_2.pdf")
    ps = []; ps2=[];
    for (i, online) in enumerate(onlines)
        t = online[:,1]./60
        p2 = scatter(t, online[:,2], label="", title="GalOx $(i)",ylabel = "AEW in nm", markersize = 3,
        markerstrokewidth = 0.7)
        vline!([t_cop[i]/60], label="copper addition", color="red")
        push!(ps2, p2)

        p = scatter(t,online.diff_AEW_raw, markersize = 3, markeralpha=0.4,
            markerstrokewidth = 0.2, label="Finite Difference",
            xlabel="Time in hours", ylim=(0,10), ylabel = "dAEW in nm/hour")
        scatter!(t,online.diff_AEW_sgol, markersize = 3, markerstrokewidth = 0.2, markershape = :diamond, label="Savitzky-Golay")
        scatter!(t,online.diff_AEW_loess, markersize = 3, markerstrokewidth = 0.2, markershape = :xcross, label="Loess")
        vline!([t_cop[i]/60], label="copper addition", color="red")
        push!(ps, p)
    end
    ps_tot = plot(vcat([p for p in ps2], [p for p in ps])..., size=(1000,550), layout=(2,Int64(length(onlines))),
        legendfontsize = 7,
        titlelocation = :left,
        bottom_margin=10Plots.px,
        left_margin=10Plots.px,
        tickfontsize = 10,
        xlabelfontsize = 10,
        ylabelfontsize = 10,
        grid = false,
        framestyle = :box,
        legend=[true false false true false false],
        ylabel = ["AEW in nm" "" "" L"$-\frac{d}{dt}$(AEW) in nm/hour" "" ""],
        xlabel = ["" "" "" "Time in hours" "Time in hours" "Time in hours"],
        title = ["GalOx 2 (0.12M GuHCl)" "GalOx 3 (0.20M GuHCl)" "GalOx 5 (0.60M GuHCl)" "" "" ""]
        )
        
    display(ps_tot)
    if save
        savefig(ps_tot, filename)
    end
    return ps_tot
end

function plot_specific_k(LDH_online, GalOx_online, tcop_adapted, GuHCl; save=false, filename="figs/specific_k.pdf")
    mshapes = [:circle, :utriangle, :star5, :diamond, :hexagon, :dtriangle, 
    :rtriangle, :ltriangle, :pentagon, :heptagon, :octagon, :star4, :star6, :star7, :star8, :vline, :hline, :+, :x]
    # pl = plot()
    # for (i,online) in enumerate(LDH_online)
    #     online.kI_meas = -online.dIdt ./ (online.I_ss.-minimum(online.I_ss))
    #     scatter!(online[:,1]./60, online.kI_meas, label="LDH $(i+3)", 
    #     xlabel="Time in hours", ylabel="kI",
    #     ylim=(0,4), legend = :bottomright,
    #     markershape = mshapes[i+3],
    #     )
    # end
    
    pg = plot()
    for (i,online) in enumerate(GalOx_online)
        online.kI_meas = -online.dIdt ./ (online.I_ss.-minimum(online.I_ss))
        scatter!((online[:,1]./60)[Float64.(online[:,1]) .< tcop_adapted[i]], 
        online.kI_meas[Float64.(online[:,1]) .< tcop_adapted[i]], 
        label="GalOx $i", 
        xlabel="Time in hours", ylabel="kI",
        ylim=(0,2), legend = :topright,
        markershape = mshapes[Int(ceil(i/2))],
        )
    end
    
    pg2 = plot()
    for (i,online) in enumerate(GalOx_online)
        scatter!((online[:,1]./60)[Float64.(online[:,1]) .> tcop_adapted[i]] .- (online[:,1]./60)[Float64.(online[:,1]) .> tcop_adapted[i]][1], 
        online.kI_meas[Float64.(online[:,1]) .> tcop_adapted[i]], 
        label="GalOx $i", 
        xlabel="Time in hours", ylabel="kI", legend = :topright,
        ylim=(0,25), 
        markershape = mshapes[Int(ceil(i/2))],
        )
    end

    pl = plot()
    onlines = []
    for (i,online) in enumerate(GalOx_online)
        online.GuHCl = [GuHCl[i] for j in 1:length(online.kI_meas)]
        #scatter!([GuHCl[i] for j in 1:length(online.kI_meas[Float64.(online[:,1]) .> tcop_adapted[i]])],
        #online.kI_meas[Float64.(online[:,1]) .> tcop_adapted[i]], 
        #label="GalOx $i", 
        #xlabel="Guanidine HCl", ylabel="kI", legend = :topright,
        #ylim=(0,50), 
        #markershape = mshapes[Int(ceil(i/2))],
        #)
        online = online[Float64.(online[:,1]) .> tcop_adapted[i],:]
        if GuHCl[i] == 0.12
            online = online[Float64.(online[:,1]) .< tcop_adapted[i].+0.2*60,:]
        elseif GuHCl[i] == 0.2
            online = online[Float64.(online[:,1]) .< tcop_adapted[i].+0.5*60,:]
        else
            online = online[Float64.(online[:,1]) .< tcop_adapted[i].+60,:]
        end
        push!(onlines, online)
    end
    df = vcat(onlines...)
    replace!(df.kI_meas, missing=>0, Inf=>0, NaN=>0)
    df = df[df.kI_meas .> 0,:]
    p = violin(string.(df.GuHCl), df.kI_meas, line = 0, fill = (0.2, :blue))
    boxplot!(string.(df.GuHCl), df.kI_meas, line = (2, :black), fill = (0.3, :orange), ylim=(0,25))
    #GalOx_online[1].dIdt[Float64.(GalOx_online[1][:,1]) .> GalOx_offline.t_cop[1]]
    pt2 = plot(pg,pg2, p, layout=(1,3), size=(1000,350),
        title=["(A) LDH" "(B) GalOx before copper" "(C) GalOx after copper"],
        ylabel = [L"$k_I$ in 1/hour" "" ""],
        xlabel = ["Time in hours" "Time in hours"],
        legendfontsize = 10,
        titlelocation = :left,
        bottom_margin=20Plots.px,
        left_margin=20Plots.px,
        tickfontsize = 10,
        xlabelfontsize = 10,
        ylabelfontsize = 10,
        grid = false,
        framestyle = :box,
        )
    display(pt2)
    if save
        savefig(pt2, filename)
    end
    return pt2
end


function plot_specific_k_violin(LDH_online, GalOx_online, tcop_adapted, GuHCl; save=false, filename="figs/specific_k.pdf")
    onlines = []
    for (i,online) in enumerate(GalOx_online)
        online.GuHCl = [GuHCl[i] for j in 1:length(online.kI_meas)]
        online = online[Float64.(online[:,1]) .< tcop_adapted[i],:]
        if GuHCl[i] == 0.12
            online = online[Float64.(online[:,1]) .< 20,:]
        elseif GuHCl[i] == 0.2
            online = online[Float64.(online[:,1]) .< 20,:]
        else
            #online = online[Float64.(online[:,1]) .< tcop_adapted[i],:]
        end
        push!(onlines, online)
    end
    df = vcat(onlines...)
    replace!(df.kI_meas, missing=>0, Inf=>0, NaN=>0)
    df = df[df.kI_meas .> 0,:]
    p = violin(string.(df.GuHCl), df.kI_meas, line = 0, fill = (0.2, :blue), label = "")
    boxplot!(string.(df.GuHCl), df.kI_meas, line = (2, :black), 
        fill = (0.3, :orange), ylim=(0,1.75), label = "")
    
    onlines = []
    for (i,online) in enumerate(GalOx_online)
        online.GuHCl = [GuHCl[i] for j in 1:length(online.kI_meas)]
        online = online[Float64.(online[:,1]) .> tcop_adapted[i],:]
        if GuHCl[i] == 0.12
            online = online[Float64.(online[:,1]) .< tcop_adapted[i].+0.2*60,:]
        elseif GuHCl[i] == 0.2
            online = online[Float64.(online[:,1]) .< tcop_adapted[i].+0.5*60,:]
        else
            online = online[Float64.(online[:,1]) .< tcop_adapted[i].+60,:]
        end
        push!(onlines, online)
    end
    df = vcat(onlines...)
    replace!(df.kI_meas, missing=>0, Inf=>0, NaN=>0)
    df = df[df.kI_meas .> 0,:]
    p2 = violin(string.(df.GuHCl), df.kI_meas, line = 0, fill = (0.2, :blue), label = "")
    boxplot!(string.(df.GuHCl), df.kI_meas, line = (2, :black), 
        fill = (0.3, :orange), ylim=(0,25), label = "", legend = :topright)
    
    pt2 = plot(p, p2, layout=(1,2), size=(650,350),
        title=["(A) GalOx before copper" "(B) GalOx after copper"],
        ylabel = [L"$k_I$ in 1/hour" L"$k_I$ in 1/hour"],
        xlabel = ["GuHCl in M" "GuHCl in M"],
        legendfontsize = 10,
        titlelocation = :left,
        bottom_margin=20Plots.px,
        left_margin=20Plots.px,
        tickfontsize = 10,
        xlabelfontsize = 10,
        ylabelfontsize = 10,
        grid = false,
        framestyle = :box,
        )
    display(pt2)
    if save
        savefig(pt2, filename)
    end
    return pt2
end