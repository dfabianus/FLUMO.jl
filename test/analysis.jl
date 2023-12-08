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

# Open-loop simulation of the models 
ps = simulate_LDH_experiment.(LDH_online, LDH_offline)

# Simulation of the direct soft sensor (method 1)
pt = simulate_LDH_soft_sensor.(LDH_online[[40,41,42]], LDH_offline[[40,41,42]])
pt2 = plot(pt..., layout=(1,3), size=(1000,550),
    title=["Process 4" "" "Process 5" "" "Process 6" ""],
    ylabel=["dAEW in nm/hour" "Concentration in g/L" "" "" "" ""],
    xlabel=["" "Time in hours" "" "Time in hours" "" "Time in hours"],
    legend=[true :bottomleft false false false false])
savefig(pt2, "figs/LDH_simulations_2.pdf")
# Simulation of the indirect UKF-based soft sensor (method 2)

plot_AEW_vs_dAEW(LDH_online[[11,12,8]], save=true, filename="figs/plot_AEW_vs_dAEW.pdf")

# save the preprocessed data
for (i,online) in enumerate(LDH_online)
    CSV.write("data/LDH/LDH_preprocessed_$i.csv", online)
end

#savefig(ptot_LDH, "figs/LDH_simulations.png")

#### -------------- GalOx ---------------------------------------------------#
diff_AEW_LDH!.(GalOx_online)

plot_multiple_AEW(GalOx_online, save=false)           # first the raw AEW
plot_multiple_dAEW_adapted(GalOx_online, save=false)  # subtract the up pulses
plot_multiple_diff_aew(GalOx_online, save=false)      # then differentiate
plot_multiple_integrals_AEW(GalOx_online, save=false) 
plot_AEW_vs_dAEW_c(GalOx_online[[1,2,4,5,6,7]], GalOx_offline[[1,2,4,5,6,7],"t_cop"], save=true)

galox_offline_r = [GalOx_offline[i,:] for i in 1:size(GalOx_offline)[1]]
pt = simulate_GalOx_soft_sensor.(GalOx_online[[1,2,4,5,6,7]], galox_offline_r[[1,2,4,5,6,7]])
pt2 = plot(pt..., layout=(2,3), size=(1000,550),
    title=["GalOx 1" "GalOx 2" "GalOx 3" "GalOx 4" "GalOx 5" "GalOx 6"],
    ylabel = ["Concentration in g/L" "" "" "Concentration in g/L" "" ""],
    xlabel = ["" "" "" "Time in hours" "Time in hours" "Time in hours"],
    legend=[false false false :topleft false false])
savefig(pt2, "figs/GalOx_simulations_2.pdf")
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