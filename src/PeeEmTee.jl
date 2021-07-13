module PeeEmTee

using HDF5
using DocStringExtensions
using StatsBase
using Distributions
using Optim
using Plots
using RecipesBase
import Base
export high_voltages, WaveSet, calculate_charges, ChargeDist, pmtresp,
       pre_fit, gauss, pmtresp_fit, pmtresp_uap, plot, plot!, calculate_transit_times, bin_data

include("tools.jl")
include("fit.jl")
include("plot.jl")

end # module
