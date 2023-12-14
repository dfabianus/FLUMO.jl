module FLUMO

export get_LDH_data, get_GalOx_data, get_HRP_data
export p_LDH
export simulate_LDH_experiment, simulate_GalOx_experiment
export simulate_LDH_soft_sensor, simulate_GalOx_soft_sensor
export corr_I0, corr_aew
export moving_average, hampel
export diff_AEW_LDH!
export plot_multiple_AEW, plot_multiple_dAEW_adapted, plot_multiple_diff_aew, plot_multiple_integrals_AEW
export plot_AEW_vs_dAEW, plot_AEW_vs_dAEW_b, plot_AEW_vs_dAEW_c, plot_AEW_vs_dAEW_d
export plot_intensity, plot_intensity_2
export plot_specific_k, plot_specific_k_violin

using XLSX
using DataFrames
#using Interpolations
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
using DataInterpolations
using SavitzkyGolay
using Loess
using StatsPlots
using BioprocessingModelLibrary.Refolding

include("utils.jl")
include("preprocessing.jl")
include("correlations.jl")
include("flumo-data-import.jl")
include("flumo_archive.jl")
include("parameters.jl")
include("model_simulations.jl")
include("plot_functions.jl")


end
