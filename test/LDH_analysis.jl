using Revise
using FLUMO
using Experiments
using BioprocessingModelLibrary.Refolding

using XLSX
using DataFrames
using Interpolations
using Random
using Distributions
using LsqFit
using Plots
using LaTeXStrings
using ModelingToolkit
using MonteCarloMeasurements
using StatsBase
using LowLevelParticleFilters, LinearAlgebra, StaticArrays, Distributions
using JLD2
using OrdinaryDiffEq

LDH_online, LDH_offline = get_LDH_data()

offline = LDH_offline[1]
online = LDH_online[1]

mP_pulses = vcat(offline.c_p[1], diff(offline.c_p[1:end-1])).*offline.V[1:end-1]./1000
mD_pulses = offline.c_GuHCl[1:end-1].*offline.V[1:end-1]./1000
V_pulses = offline.V[1:end-1]./1000
pulse_times = offline.time[1:end-1]./60
end_time = online[end,1]./60

# pulse_times = [1.0, 2.0, 3.0, 4.0, 5.0]
# mP_pulses = [0.1, 0.2, 0.2, 0.5, 0.1]
# mD_pulses = [0.1, 0.1, 0.1, 0.1, 0.1]
# mC_pulses = [0.4, 0.2, 0.2, 0.1, 0.3]
# V_pulses = [0.2, 0.2, 0.2, 0.2, 0.2]
aₙ = 1.3343 ± 1.5837
aₐ = 12.0465 ± 1.3641
bₙ = -8.6824 ± 0.7182
bₐ = -16.7869 ± 2.5716
sys = FLUMO_COMBINED_6(pulse_times[2:end], mP_pulses[2:end], mD_pulses[2:end], zeros(5), V_pulses[2:end])
p = (sys.reactor.c_Din => 5, sys.kinetics.a_n => aₙ, sys.kinetics.b_n => bₙ, 
    sys.kinetics.a_a => aₐ, sys.kinetics.b_a => bₐ, sys.kinetics.a_ac => 1.0, sys.kinetics.a_ic => 1.0,
    sys.kinetics.a_nc => 1.0, sys.kinetics.a_cn => 1.0
)
u0 = [sys.reactor.D => mD_pulses[1], sys.reactions_ODE.I => mP_pulses[1], sys.reactions_ODE.A => 0.0, 
    sys.reactions_ODE.C => 0.0, sys.reactions_ODE.IC => 0.0, sys.reactions_ODE.NC => 0.0,
    sys.reactions_ODE.N => 0.0, sys.reactor.V => V_pulses[1].±0.001
]
tspan = (0.,5.)
oprob = ODEProblem(sys, u0, tspan, p)
osol  = solve(oprob, Tsit5())
ts = range(tspan..., length=300)

p = plot(ts, osol(ts, idxs=sys.reactions_ODE.N).u, label = "N(t)", 
title = "Testset MTK Catalyst Connections", xlabel="Time (h)", ylabel="Concentration (mol/L)",
ylims=(0,0.5))
p = plot!(ts, osol(ts, idxs=sys.reactions_ODE.I).u, label = "I(t)")
p = plot!(ts, osol(ts, idxs=sys.reactions_ODE.A).u, label = "A(t)")
scatter!(offline.time./60, offline[!,"c_N "], label = "N(t) data")
scatter!(offline.time./60, offline.c_p, label = "P(t) data")
p = plot(ts, osol(ts, idxs=sys.kinetics.k_n).u, label = "k_n", title = "Testset MTK Catalyst Connections", xlabel="Time (h)", ylabel="Concentration (mol/L)")