module FLUMO

export get_LDH_data

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

using Experiments
using BioprocessingModelLibrary

include("flumo-data-import.jl")
include("flumo1.jl")

end
