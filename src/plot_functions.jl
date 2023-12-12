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
    #annotate!(1.75, 2.7E5, text("Insoluble \nAggregates", 4, :right, 15))
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
        push!(RMSE_raw, sqrt(mean((online.dAEW_adapted[2:end] .- (-cumsum(online.diff_AEW_raw[2:end].*diff(online[:,1]./60)))).^2)))
        push!(RMSE_sgol, sqrt(mean((online.dAEW_adapted[2:end] .- (-cumsum(online.diff_AEW_sgol[2:end].*diff(online[:,1]./60)))).^2)))
        push!(RMSE_loess, sqrt(mean((online.dAEW_adapted[2:end] .- (-cumsum(online.diff_AEW_loess[2:end].*diff(online[:,1]./60)))).^2)))
    end
    ps_tot = plot([p for p in ps]...,size=(3000,3000))
    display(ps_tot)
    println("sum RMSE raw: ", sum(RMSE_raw))
    println("sum RMSE sgol: ", sum(RMSE_sgol))
    println("sum RMSE loess: ", sum(RMSE_loess))
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
            ylabel!(L"$\frac{d}{dt}$(AEW) in nm/hour")
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
        legend=[:topleft false false :topleft false false],
        ylabel = ["AEW in nm" "" "" L"$\frac{d}{dt}$(AEW) in nm/hour" "" ""],
        xlabel = ["" "" "" "Time in hours" "Time in hours" "Time in hours"],
        )
    display(ps_tot)
    if save
        savefig(ps_tot, filename)
    end
    return ps_tot
end