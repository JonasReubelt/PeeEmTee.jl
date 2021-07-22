module PeeEmTee

using HDF5
using DocStringExtensions
using StatsBase
using Distributions
using Optim
using RecipesBase
using DelimitedFiles
import Base
export high_voltages, WaveSet, calculate_charges, ChargeDist, pmtresp,
       pre_fit, gauss, norm_gauss, pmtresp_fit, pmtresp_uap, plot, plot!,
       calculate_transit_times, bin_data, subtract_baseline, simulate_charges,
       baseline_max, PMTRespFit, PMTRespUapFit, make_qfunc, PreFit, mean_waveform

include("tools.jl")
include("fit.jl")
include("plot.jl")

end # module
