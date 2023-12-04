using Revise
using FLUMO
using BioprocessingModelLibrary.Refolding

using OrdinaryDiffEq
using ModelingToolkit

# using XLSX
# using DataFrames
# using Interpolations
# using Random
# using Distributions
# using LsqFit
using Plots
# using LaTeXStrings
# using ModelingToolkit
using MonteCarloMeasurements
# using StatsBase
# using LowLevelParticleFilters, LinearAlgebra, StaticArrays, Distributions
# using JLD2
# using OrdinaryDiffEq
using Loess
using SavitzkyGolay
#---------------------------------------------------------------------------#

LDH_online, LDH_offline = get_LDH_data()
GalOx_online, GalOx_offline = get_GalOx_data()
HRP_online = get_HRP_data()

#### --------------- derivatives  ------------------------------------------#
scatter(LDH_online[1][2:end,1]./60, LDH_online[1].Integral)
scatter(LDH_online[1][2:end,1]./60, LDH_online[1][!,2])

function diff_AEW_LDH!(online)
    t = deepcopy(online[:,1]./60)
    #F = deepcopy(online.Integral)
    AEW = deepcopy(Float64.(online[:,2]))
    dAEWdt = diff(AEW)./diff(t)
    idxInf = findall(x->x>5, dAEWdt)
    for i in idxInf
        AEW[i+1:end] = AEW[i+1:end].-diff(AEW)[i]
    end
    #AEW = hampel(AEW, 1, 3)

    dAEWdt = diff(AEW)./diff(t)
    dAEWdt_smooth = savitzky_golay(AEW, 15, 5, deriv=1).y[2:end]./diff(t)
    for i in idxInf
        dAEWdt[i] = (dAEWdt[i+1] + dAEWdt[i-1])/2
        dAEWdt_smooth[i] = (dAEWdt_smooth[i+1] + dAEWdt_smooth[i-1])/2
    end
    dAEWdt = hampel(dAEWdt, 2, 2)
    dAEWdt_smooth = hampel(dAEWdt_smooth, 2, 2)
    dAEWdt_smooth_loess = predict(loess(t[2:end], dAEWdt, span=0.09), t[2:end])
    dAEWdt = [dAEWdt[1]; dAEWdt]
    dAEWdt_smooth = [dAEWdt_smooth[1]; dAEWdt_smooth]
    dAEWdt_smooth_loess = [dAEWdt_smooth_loess[1]; dAEWdt_smooth_loess]

    p1 = scatter(t,AEW.-AEW[1])
    p = scatter(t,dAEWdt)
    plot!(t,dAEWdt_smooth, linewidth=4)
    plot!(t,dAEWdt_smooth_loess)
    online.diff_AEW = [maximum([0, -d]) for d in dAEWdt_smooth_loess]
    return p1, p
end

ps3 = diff_AEW_LDH!.(LDH_online)
ints = -[cumsum(online.diff_AEW[2:end].*diff(online[:,1]./60)) for online in LDH_online]
plot([p[1] for p in ps3]...,size=(3000,3000))
plot([scatter(online.diff_AEW) for online in LDH_online]...,size=(3000,3000))
plot([scatter(online[:,1]./60, online[!,2]) for online in LDH_online]...,size=(3000,3000))
for (i,online) in enumerate(LDH_online)
    p = scatter([p[1] for p in ps3][i])
    scatter!(online[2:end,1]./60, ints[i])
    display(p)
end

#### --------------- LDH ---------------------------------------------------#

# output correlations to fluorescence intensity and average emission wavelength
fit_u1, confInt, p_u1 = corr_I0(LDH_offline[1:38])
fit_u2, p_u2 = corr_aew(LDH_offline[1:38])

ps = simulate_LDH_experiment.(LDH_online, LDH_offline)
@register_symbolic F_fun(t)
@register_symbolic dAEWdt_fun(t)
pt = simulate_LDH_soft_sensor.(LDH_online, LDH_offline)

ptot_LDH = plot([p[2] for p in ps]..., layout=(7,7), size=(3000,3000), legendfontsize = 10,
        title = "",
        titlelocation = :left,
        bottom_margin=30Plots.px,
        left_margin=20Plots.px,
        tickfontsize = 10,
        xlabelfontsize = 10,
        ylabelfontsize = 10,
        grid = false,
        framestyle = :box)

savefig(ptot_LDH, "figs/LDH_simulations.png")

#### -------------- GalOx ---------------------------------------------------#
p = []
for (idx,online) in enumerate(GalOx_online)
    t = online[:,"Process Time [min]"]./60
    I = online.Integral
    push!(p, plot(t,I))
end
plot(p...)

p = []
for (idx,online) in enumerate(GalOx_online)
    t = online[:,"Process Time [min]"]./60
    AEW = online[:,"Average emission wavelength [nm]"]
    push!(p, plot(t,AEW))
end
plot(p...)

galox_offline_r = [GalOx_offline[i,:] for i in 1:size(GalOx_offline)[1]]
ps2 = simulate_GalOx_experiment.(GalOx_online, galox_offline_r)

ptot_GalOx = plot([p[2] for p in ps2]..., layout=(3,4), size=(2000,1500), legendfontsize = 10,
        title = "",
        titlelocation = :left,
        bottom_margin=30Plots.px,
        left_margin=20Plots.px,
        tickfontsize = 10,
        xlabelfontsize = 10,
        ylabelfontsize = 10,
        grid = false,
        framestyle = :box)

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