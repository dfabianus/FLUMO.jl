using Revise
using FLUMO
using BioprocessingModelLibrary.Refolding

using OrdinaryDiffEq
using ModelingToolkit

using Plots
using MonteCarloMeasurements
using CSV

#---------------------------------------------------------------------------#

LDH_online, LDH_offline = get_LDH_data()
GalOx_online, GalOx_offline = get_GalOx_data()
HRP_online = get_HRP_data()

#### --------------- LDH ---------------------------------------------------#

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
ps = simulate_LDH_experiment.(LDH_online, LDH_offline)

# Simulation of the direct soft sensor (method 1)
exps = [40,41,42]
exps = [11,12,8]
exps = [11,12,8,40,41,42]
exps = [14 16 32 33 34 37]
exps = [13,14,16,40,41,42]
exps = collect(1:42)

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
diff_AEW_LDH!.(GalOx_online)

plot_multiple_AEW(GalOx_online, save=false)           # first the raw AEW
plot_multiple_dAEW_adapted(GalOx_online, save=false)  # subtract the up pulses
plot_multiple_diff_aew(GalOx_online, save=false)      # then differentiate
plot_multiple_integrals_AEW(GalOx_online, save=false) 
plot_AEW_vs_dAEW_c(GalOx_online[[1,2,4,5,6,7]], GalOx_offline[[1,2,4,5,6,7],"t_cop"], save=true)

# Soft-sensor validation 
galox_offline_r = [GalOx_offline[i,:] for i in 1:size(GalOx_offline)[1]]
pt = simulate_GalOx_soft_sensor.(GalOx_online[[1,2,4,5,6,7]], galox_offline_r[[1,2,4,5,6,7]]; distinct=false)
pt2 = plot(pt..., layout=(2,3), size=(1000,550),
    title=["GalOx 1" "GalOx 2" "GalOx 3" "GalOx 4" "GalOx 5" "GalOx 6"],
    ylabel = ["Concentration in g/L" "" "" "Concentration in g/L" "" ""],
    xlabel = ["" "" "" "Time in hours" "Time in hours" "Time in hours"],
    legend=[false false false :topleft false false])
savefig(pt2, "figs/GalOx_val_1.pdf")

# Soft-sensor validation with distinct N and A --> BAD
pt = simulate_GalOx_soft_sensor.(GalOx_online[[1,2,4,5,6,7]], galox_offline_r[[1,2,4,5,6,7]]; distinct=true)
pt2 = plot(pt..., layout=(2,3), size=(1000,550),
    title=["GalOx 1" "GalOx 2" "GalOx 3" "GalOx 4" "GalOx 5" "GalOx 6"],
    ylabel = ["Concentration in g/L" "" "" "Concentration in g/L" "" ""],
    xlabel = ["" "" "" "Time in hours" "Time in hours" "Time in hours"],
    legend=[false false false :topleft false false])
savefig(pt2, "figs/GalOx_val_2.pdf")







#### --------------- HRP ---------------------------------------------------#

p = []
for (idx,online) in enumerate(HRP_online)
    t = online[:,"Process Time [h]"]./60
    I = online[!, 3]
    push!(p, plot(t,I))
end
plot(p...)

p = []
for (idx,online) in enumerate(HRP_online)
    t = online[:,"Process Time [h]"]./60
    AEW = online[!,2]
    push!(p, plot(t,AEW))
end
plot(p...)