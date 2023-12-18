using Revise
using FLUMO
using BioprocessingModelLibrary.Refolding

using OrdinaryDiffEq
using ModelingToolkit

using Plots
using MonteCarloMeasurements
using CSV
using LaTeXStrings

#---------------------------------------------------------------------------#

LDH_online, LDH_offline = get_LDH_data()
GalOx_online, GalOx_offline = get_GalOx_data()
HRP_online = get_HRP_data()

#### --------------- LDH ---------------------------------------------------#

# correlation AEW to total reaction rate 
CSV.read("data/dAEW_NA_data.csv", DataFrame)

# numeric differentiation of AEW
diff_AEW_LDH!.(LDH_online) 

# output correlations to fluorescence intensity and average emission wavelength
fit_u1, confInt, p_u1 = corr_I0(LDH_offline[1:38])
fit_u2, p_u2 = corr_aew(LDH_offline[1:38])

# plot the AEWs and their derivatives
plot_multiple_AEW(LDH_online, save=false)           # first the raw AEW
plot_multiple_dAEW_adapted(LDH_online, save=false)  # subtract the up pulses
plot_multiple_diff_aew(LDH_online, save=false)      # then differentiate

# When reintegrate the obtained derivative, we see different degrees of errors 
plot_multiple_integrals_AEW(LDH_online, save=false) 

# Illustration of derivative 
plot_AEW_vs_dAEW(LDH_online[[11,12,8]], save=true, filename="figs/plot_AEW_vs_dAEW.pdf")

# Open-loop simulation of the models 
ps = simulate_LDH_experiment.(LDH_online[[1,3,11,21]], LDH_offline[[1,3,11,21]])
ps2 = plot([p[2] for p in ps]..., size=(1200,350),
    layout=(1,4),
    ylabel=["Intensity" "" "" ""],
    xlabel="Time [h]",
    title=["(C)" "(D)" "(E)" "(F)"],
    legendfontsize = 10,
    titlelocation = :left,
    bottom_margin=20Plots.px,
    left_margin=20Plots.px,
    tickfontsize = 10,
    xlabelfontsize = 10,
    ylabelfontsize = 10,
    grid = false,
    framestyle = :box,
    legend = :bottomright)
savefig(ps2, "figs/FluMo1_openloop_intensity.pdf")

# Soft-sensor validation 
exps = [13,14,16,40,41,42]
pt = simulate_LDH_soft_sensor.(LDH_online[exps], LDH_offline[exps],distinct=false)
pt2 = plot(pt..., layout=(2,3), size=(1000,550),
    title=["LDH 4" "LDH 5" "LDH 6" "LDH 7" "LDH 8" "LDH 9"],
    ylabel = ["Concentration in g/L" "" "" "Concentration in g/L" "" ""],
    xlabel = ["" "" "" "Time in hours" "Time in hours" "Time in hours"],
    xlim = (0, 2.6),
    legend=[false false false :topleft false false])
savefig(pt2, "figs/LDH_val_1.pdf")

# Soft-sensor validation with distinct N and A
pt = simulate_LDH_soft_sensor.(LDH_online[exps], LDH_offline[exps],distinct=true)
pt2 = plot(pt..., layout=(2,3), size=(1000,550),
    title=["LDH 4" "LDH 5" "LDH 6" "LDH 7" "LDH 8" "LDH 9"],
    ylabel = ["Concentration in g/L" "" "" "Concentration in g/L" "" ""],
    xlabel = ["" "" "" "Time in hours" "Time in hours" "Time in hours"],
    xlim = (0, 2.6),
    legend=[false false false :topleft false false])
savefig(pt2, "figs/LDH_val_2.pdf")

#### -------------- GalOx ---------------------------------------------------#
exps_galox = [3,6,1,7,2,5]
diff_AEW_LDH!.(GalOx_online)

plot_multiple_AEW(GalOx_online, save=false)           # first the raw AEW
plot_multiple_dAEW_adapted(GalOx_online, save=false)  # subtract the up pulses
plot_multiple_diff_aew(GalOx_online, save=false)      # then differentiate
plot_multiple_integrals_AEW(GalOx_online, save=false) 
plot_AEW_vs_dAEW_c(GalOx_online[exps_galox], GalOx_offline[exps_galox,"t_cop"], save=false)
plot_AEW_vs_dAEW_d(GalOx_online[[6,1,2]], GalOx_offline[[6,1,2],"t_cop"], save=true)
plot_intensity(LDH_online[41]; save=true, filename = filename="figs/intensity_illustration.pdf")

# Soft-sensor validation 
galox_offline_r = [GalOx_offline[i,:] for i in 1:size(GalOx_offline)[1]]
pt = simulate_GalOx_soft_sensor.(GalOx_online[exps_galox], galox_offline_r[exps_galox]; distinct=false)
pt2 = plot(pt..., layout=(2,3), size=(1000,550),
    title  = ["GalOx 1" "GalOx 2" "GalOx 3" "GalOx 4" "GalOx 5" "GalOx 6"],
    ylabel = ["Concentration in g/L" "" "" "Concentration in g/L" "" ""],
    xlabel = ["" "" "" "Time in hours" "Time in hours" "Time in hours"],
    legend=[false false false false false :left])
savefig(pt2, "figs/GalOx_val_1.pdf")

# Soft-sensor validation with distinct N and A --> BAD
pt = simulate_GalOx_soft_sensor.(GalOx_online[exps_galox], galox_offline_r[exps_galox]; distinct=true)
pt3a = plot_intensity_2(GalOx_online[4]; save=false, filename = filename="figs/intensity_illustration.pdf")
pt3b = plot(pt[3], title = "(B) State estimation of GalOx 3", ylabel = "Concentration in g/L", xlabel = "Time in hours")
pt3 = plot(pt3a, pt3b, layout=(1,2), size=(1000,350),
bottom_margin=20Plots.px,
left_margin=20Plots.px)
savefig(pt3, "figs/GalOx_val_2.pdf")

#### --------------- Mechanistic modeling ------------------------------------#
pm = simulate_LDH_soft_sensor!.(LDH_online[[12]], LDH_offline[[12]],distinct=true)
#pm2 = plot_AEW_vs_dAEW(LDH_online[[12]], save=false)
pt2 = plot([pm[1], pm2]..., layout=(2,3), size=(1000,550),
    title=["LDH 4" "LDH 5" "LDH 6" "LDH 7" "LDH 8" "LDH 9"],
    ylabel = ["Concentration in g/L" "" "" "Concentration in g/L" "" ""],
    xlabel = ["" "" "" "Time in hours" "Time in hours" "Time in hours"],
    xlim = (0, 2.6),
    legend=[false false false :topleft false false])

# specific kI over time LDH
plot_specific_k_violin(LDH_online[exps_galox], GalOx_online[exps_galox],GalOx_offline.t_cop[exps_galox], 
    GalOx_offline.GuHCl[exps_galox]; save=true)