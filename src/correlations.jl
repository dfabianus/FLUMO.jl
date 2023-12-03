function corr_I0(offline)
    cp = []; I0 = [];
    for df_offline in offline
        append!(I0, df_offline.I0)
        append!(cp, df_offline.c_p)
    end
    monod = (x,p) -> p[1] .* x ./ (p[2] .+ x) # monod
    monod_fit = LsqFit.curve_fit(monod, cp, I0, [100000.0, 100000.0])
    mpred = [monod(x,monod_fit.param) for x in [0:0.1:1]]
    confidence_intervals_68 = confidence_interval(monod_fit, 0.32)
    # Plot the correlation of I0 and cp
    lower_bounds_68 = [monod(x,[confidence_intervals_68[1][1], confidence_intervals_68[2][2]]) for x in 0:0.1:1]
    upper_bounds_68 = [monod(x,[confidence_intervals_68[1][2], confidence_intervals_68[2][1]]) for x in 0:0.1:1]
    p=plot([0:0.1:1], lower_bounds_68, fillranges = upper_bounds_68, fillalpha = 0.25, 
        label = L"Confidence band ($1\cdot\sigma$)", fillcolor = 2, linewidth=0)
    plot!([0:0.1:1], mpred, color = 2, linewidth = 2, label = "model fit")
    scatter!(cp, I0, xlabel = L"$c_p$ [g/L]", ylabel = L"$I_0$ [a.u]", 
        color = 1, markersize = 3, label = "Training data",
        legendfontsize = 10,
        tickfontsize = 10,
        xlabelfontsize = 10,
        ylabelfontsize = 10,
        grid = true,
        framestyle = :box)
    display(p)
    return (fit=monod_fit, confInt=confidence_intervals_68, plot=p)
end

function corr_aew(dfs_offline)
    ### cN vs delta aew
    daews = []
    cns = []
    cp = []
    ca = []
    c_GuHCl = []
    for i in range(13,38)
        if dfs_offline[i][!,end][end] + dfs_offline[i].c_A[end] > 0.0  #&& dfs_offline[i].c_p[end] > 0.2 && dfs_offline[i][!,end][end] + dfs_offline[i].c_A[end] < 0.5 #&& dfs_offline[i].delta_aew[end] > 0.2 && dfs_offline[i].delta_aew[end] < 0.95
            push!(daews, dfs_offline[i].delta_aew[end])
            push!(cns, dfs_offline[i][!,end][end])
            push!(cp, dfs_offline[i].c_p[end])
            push!(c_GuHCl, dfs_offline[i].c_GuHCl[end])
            push!(ca, dfs_offline[i].c_A[end])
        end
    end
    #model_data = DataFrame(cnas = Float64.(ca.+cns), ca = Float64.(ca), cns = Float64.(cns), daews = Float64.(daews), cp = Float64.(cp), c_GuHCl = Float64.(c_GuHCl))
    p=scatter((cns.+ca)./cp, daews, xlabel = "(cn+ca)/cp", ylabel = "daew", label = "data", 
        title = L"$\Delta$aew vs. spec. native protein concentration",
        legendfontsize = 10,
        tickfontsize = 10,
        xlabelfontsize = 10,
        ylabelfontsize = 10,
        grid = true,
        framestyle = :box)
    monod = (x,p) -> abs(p[1]) .* x ./ (abs(p[2]) .+ x) .+ abs(p[3]).*x # monod
    monod_std = (x,p) -> (p[1]) .* x ./ ((p[2]) .+ x) # monod
    monod_fit = LsqFit.curve_fit(monod, (cns.+ca)./cp, daews, [1.0, 0.5, 0.5])
    monod_fit_std = LsqFit.curve_fit(monod_std, (cns.+ca)./cp, daews, [1.0, 0.5])
    monod_fit_std.param[2] = 0.2
    mpred_std = [monod_std(x,monod_fit_std.param) for x in [0:0.01:0.85]]
    monod_fit.param[2] = 0.2
    mpred2 = [monod(x,monod_fit.param) for x in [0:0.01:0.85]]
    plot!([0:0.01:0.85], mpred_std, linewidth = 2, label = "model fit")
    plot!([0:0.01:0.85], mpred2, linewidth = 2, label = "model fit 2")
    display(p)
    return (fit = monod_fit, plot=p)
end