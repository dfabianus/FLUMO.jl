# # save the preprocessed data
# for (i,online) in enumerate(LDH_online)
#     CSV.write("data/LDH/LDH_preprocessed_$i.csv", online)
# end
# ptot_LDH = plot([p[2] for p in ps]..., layout=(7,7), size=(3000,3000), legendfontsize = 10,
#         title = "",
#         titlelocation = :left,
#         bottom_margin=30Plots.px,
#         left_margin=20Plots.px,
#         tickfontsize = 10,
#         xlabelfontsize = 10,
#         ylabelfontsize = 10,
#         grid = false,
#         framestyle = :box)

function splitdf(df, pct)
    @assert 0 <= pct <= 1
    ids = collect(axes(df, 1))
    shuffle!(ids)
    sel = ids .<= nrow(df) .* pct
    return view(df, sel, :), view(df, .!sel, :)
end

function split_train_test(dfs_online, dfs_offline, experiments, train_data_percentage = 0.6)
    train_set_ids = (sample(experiments, Int(round(train_data_percentage*length(experiments))), replace = false))
    test_set_ids = setdiff(experiments, train_set_ids)
    train_set_offline = dfs_offline[train_set_ids]
    test_set_offline = dfs_offline[test_set_ids]
    train_set_online = dfs_online[train_set_ids]
    test_set_online = dfs_online[test_set_ids]
    return (
        dfs_offline_train=train_set_offline, 
        dfs_offline_test=test_set_offline, 
        dfs_online_train=train_set_online, 
        dfs_online_test=test_set_online)
end

function calculate_delta_aew!(dfs)
    for df in dfs
        aew = df[!,2]
        df.Δaew = aew .- aew[1]
        #df.ddaew = vcat([0.0], diff(df.Δaew))
        #println(df.Δaew[df.ddaew.>1.0])
        #df.Δaew[df.ddaew.>0.1] .= df.Δaew[df.ddaew.>0.1] .- df.ddaew[df.ddaew.>0.1]
       # df.ddaew_k = kalman_state_derivative(df[!,1]./60, df.ddaew, [-2.5 0.0]', [10.0 0.0; 0.0 10.0]; Q = [0.1 0.0; 0.0 0.0], R = 0.1)[1]
    end
end

function predicted_vs_observed(data, sols_train, sols_test)
    data_train = []
    data_test = []
    scatter(data.dfs_online_train[1].Integral, Measurements.value.(sols_train[1](data.dfs_online_train[1][!,1]./60, idxs = (Intensity_mon))[Intensity_mon]), 
        color=1, markersize=3, markerstrokewidth=0.5, label="training data")
    for (df,sol) in zip(data.dfs_online_train, sols_train)
        push!(data_train, DataFrame(Integral = df.Integral, Predicted = Measurements.value.(sol(df[!,1]./60, idxs = (Intensity_mon))[Intensity_mon])))
    end
    scatter!(data.dfs_online_test[1].Integral, Measurements.value.(sols_test[1](data.dfs_online_test[1][!,1]./60, idxs = (Intensity_mon))[Intensity_mon]), 
    color=2, markersize=3, markerstrokewidth=0.5, label="test data")
    for (df,sol) in zip(data.dfs_online_test, sols_test)
        push!(data_test, DataFrame(Integral = df.Integral, Predicted = Measurements.value.(sol(df[!,1]./60, idxs = (Intensity_mon))[Intensity_mon])))
    end
    data_train=Float64.(vcat(data_train...))
    data_test=Float64.(vcat(data_test...))
    p=scatter(data_train.Integral, data_train.Predicted, color=1, markersize=3, markerstrokewidth=0.5, label="training data")
    scatter!(data_test.Integral, data_test.Predicted, color=2, markersize=3, markerstrokewidth=0.5, label="test data")
    plot!([0:1000:7.5e5], x -> x, color=:black, linewidth=2, label="", 
    xlabel="observed fluorescence [a.u.]", ylabel="predicted fluorescence [a.u.]",
    title="Predicted vs. observed fluorescence")
    display(p)
    #model = lm(@formula(Predicted ~ Integral), data_train)
end

function model_I0_cp(dfs_offline_train, dfs_offline_test)
    I0_test = []; cp_test = []; I0_train = []; cp_train = []
    for df_offline_train in dfs_offline_train
        append!(I0_train, df_offline_train.I0)
        append!(cp_train, df_offline_train.c_p)
    end
    for df_offline_test in dfs_offline_test
        append!(I0_test, df_offline_test.I0)
        append!(cp_test, df_offline_test.c_p)
    end
    monod = (x,p) -> p[1] .* x ./ (p[2] .+ x) # monod
    monod_fit = LsqFit.curve_fit(monod, cp_train, I0_train, [100000.0, 100000.0])
    mpred = [monod(x,monod_fit.param) for x in [0:0.1:1]]
    confidence_intervals_68 = confidence_interval(monod_fit, 0.32)
    # Plot the correlation of I0 and cp
    lower_bounds_68 = [monod(x,[confidence_intervals_68[1][1], confidence_intervals_68[2][2]]) for x in 0:0.1:1]
    upper_bounds_68 = [monod(x,[confidence_intervals_68[1][2], confidence_intervals_68[2][1]]) for x in 0:0.1:1]
    p=plot([0:0.1:1], lower_bounds_68, fillranges = upper_bounds_68, fillalpha = 0.25, label = L"Confidence band ($1\cdot\sigma$)", fillcolor = 2, linewidth=0)
    #p=plot([0:0.1:1], lower_bounds_68, fillranges = upper_bounds_68, fillalpha = 0.25, label = L"Confidence band ($1\cdot\sigma$)", fillcolor = 2, linewidth=0)
    plot!([0:0.1:1], mpred, color = 2, linewidth = 2, label = "model fit")
    scatter!(cp_train, I0_train, xlabel = L"$c_p$ [g/L]", ylabel = L"$I_0$ [a.u]", color = 1, markersize = 3, label = "Training data")
    scatter!(cp_test, I0_test, color = 3, markersize = 3, label = "Test data",
    title = "Intensity vs. protein concentration",
    legendfontsize = 10,
    tickfontsize = 10,
    xlabelfontsize = 10,
    ylabelfontsize = 10,
    #size = (800,600),
    grid = true,
    framestyle = :box)
    display(p)
    return (fit=monod_fit, confInt=confidence_intervals_68, plot=p)
end

function model_aew_cn(dfs_offline)
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
    model_data = DataFrame(cnas = Float64.(ca.+cns), ca = Float64.(ca), cns = Float64.(cns), daews = Float64.(daews), cp = Float64.(cp), c_GuHCl = Float64.(c_GuHCl))
    
    #fm = @formula(daews ~ cnas)
    #linear_regressor = lm(fm, model_data)
    p=scatter((cns.+ca)./cp, daews, xlabel = "(cn+ca)/cp", ylabel = "daew", label = "data", 
    title = L"$\Delta$aew vs. spec. native protein concentration",
    legendfontsize = 10,
    tickfontsize = 10,
    xlabelfontsize = 10,
    ylabelfontsize = 10,
    #legend  = :bottomright,
    #size = (800,600),
    grid = true,
    framestyle = :box)
    monod = (x,p) -> abs(p[1]) .* x ./ (abs(p[2]) .+ x) .+ abs(p[3]).*x # monod
    monod_std = (x,p) -> (p[1]) .* x ./ ((p[2]) .+ x) # monod
    #monod = (x,p) -> p[1] .* x ./ (p[2] .+ x .* (1 .+ x ./ p[3]))
    #lf2 = ((cns.+ca)./cp)[((cns.+ca)./cp).<0.034761028197956777]
    #daews2 = daews[((cns.+ca)./cp).<0.034761028197956777]
    #linear_fit_2 = daews2[1]/lf2[1]
    #linear = (x,p) -> p[1] .* x .+ daews2[1]#linear
    #linear = (x,p) -> p[1] .* x .+ p[2]#linear
    # logfun = (x,p) -> p[1] .* log.(x) .+ p[2]
monod_fit = LsqFit.curve_fit(monod, (cns.+ca)./cp, daews, [1.0, 0.5, 0.5])
monod_fit_std = LsqFit.curve_fit(monod_std, (cns.+ca)./cp, daews, [1.0, 0.5])
    #lf1 = ((cns.+ca)./cp)[((cns.+ca)./cp).>0.034761028197956777]
    #daews1 = daews[((cns.+ca)./cp).>0.034761028197956777]
    #linear_fit_1 = LsqFit.curve_fit(linear, ((cns.+ca)./cp), daews, [10.0,10.0])
    
    # log_fit = LsqFit.curve_fit(logfun, cns./cp, daews, [10.0, 10.0])
    #println(monod_fit.param)
monod_fit_std.param[2] = 0.2
mpred_std = [monod_std(x,monod_fit_std.param) for x in [0:0.01:0.85]]
monod_fit.param[2] = 0.2
mpred2 = [monod(x,monod_fit.param) for x in [0:0.01:0.85]]
    #linearpred = [linear(x,linear_fit_1.param) for x in [0:0.01:0.85]]
    # logpred = [logfun(x,log_fit.param) for x in [0:0.1:0.7]]
#confidence_intervals_68 = confidence_interval(monod_fit, 0.32)
    #confidence_intervals_68_linear = confidence_interval(linear_fit, 0.32)
    # confidence_intervals_68_log = confidence_interval(log_fit, 0.32)
    # Plot the correlation of I0 and cp
#lower_bounds_68 = [monod(x,[confidence_intervals_68[1][1], confidence_intervals_68[2][2], confidence_intervals_68[3][1]]) for x in 0:0.01:0.85]
#upper_bounds_68 = [monod(x,[confidence_intervals_68[1][2], confidence_intervals_68[2][1], confidence_intervals_68[3][2]]) for x in 0:0.01:0.85]
    #lower_bounds_68_linear = [linear(x,[confidence_intervals_68_linear[1][1]]) for x in 0:0.1:0.7]
    #upper_bounds_68_linear = [linear(x,[confidence_intervals_68_linear[1][2]]) for x in 0:0.1:0.7]
    # lower_bounds_68_log = [monod(x,[confidence_intervals_68_log[1][2], confidence_intervals_68_log[2][1]]) for x in 0:0.1:0.7]
    # upper_bounds_68_log = [monod(x,[confidence_intervals_68_log[1][2], confidence_intervals_68_log[2][1]]) for x in 0:0.1:0.7]

#plot!([0:0.01:0.85], lower_bounds_68, fillranges = upper_bounds_68, fillalpha = 0.25, label = L"Confidence band ($1\cdot\sigma$)", fillcolor = 2, linewidth=0)
plot!([0:0.01:0.85], mpred_std, linewidth = 2, label = "model fit")
plot!([0:0.01:0.85], mpred2, linewidth = 2, label = "model fit 2")
    # plot!([0:0.1:0.7], lower_bounds_68_log, fillranges = upper_bounds_68_log, fillalpha = 0.25, label = L"Confidence band ($1\cdot\sigma$)", fillcolor = 3, linewidth=0)
    # plot!([0:0.1:0.7], logpred, color = 3, linewidth = 2, label = "Log fit")
    #plot!([0:0.1:0.7], lower_bounds_68_linear, fillranges = upper_bounds_68_linear, fillalpha = 0.25, label = L"Confidence band ($1\cdot\sigma$)", fillcolor = 4, linewidth=0)
    #plot!([0:0.01:0.034761028197956777], linear_fit_2.*[0:0.01:0.034761028197956777], color = 4, linewidth = 2, label = "Linear fit 1")
    #plot!([0:0.01:0.85], linearpred, color = 4, linewidth = 2, label = "Linear fit")
    display(p)
    #return (fit=monod_fit, confInt=confidence_intervals_68, plot=p, glmodel=linear_regressor)
    return (fit = monod_fit, plot=p)
end

# plotbool = true

# path = "03_Exchange_CIG_FMU/Puls1/231002_Data_Pulse1_offline.xlsx"
# dfs = []
# for i in range(1,8)
#     push!(dfs, DataFrame(XLSX.readtable(path, "run_$i", "A:J")))
# end
# df = vcat(dfs...)

# scatter(df.c_p, df.I0)

# X = Float64.(hcat(ones(length(df.c_p)), df.c_p))
# y = Float64.(df.I0)

# p_lin = X \ y

# linear = (x,p) -> p[1] .+ p[2] .* x #linear
# f2 = (x,p) -> p[1] .+ p[2] .* x .+ p[3] .* x.^2.0 # quadratic
# monod = (x,p) -> p[1] .* x ./ (p[2] .+ x) # monod

# linear_fit = LsqFit.curve_fit(linear, 
#                     df.c_p[df.c_p .< 0.5],
#                     df.I0[df.c_p .< 0.5],
#                     [100000.0, 100000.0])
# monod_fit = LsqFit.curve_fit(monod, 
#                     df.c_p,
#                     df.I0,
#                     [100000.0, 100000.0])


# if plotbool == true
#     scatter(df.c_p, df.I0, label = "Data", xlabel = L"Protein concentration $c_p$ [g/L]", ylabel = L"Intensity $I_0$", title="Intensity correlations")
#     plot!([0:0.1:1], x -> linear(x,linear_fit.param), label = "Linear fit")
#     plot!([0:0.1:1], x -> monod(x,monod_fit.param), label = "Monod fit")
# end