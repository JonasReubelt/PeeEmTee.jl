"""
$(SIGNATURES)
reads the high voltages of the wavesets saved in the hdf5 file
# Arguments
- `filename`: filename of h5 file
"""
function high_voltages(filename)
    h5 = h5open(filename)
    high_voltages = keys(h5)
    close(h5)
    high_voltages
end



struct WaveSet
    waveforms
    v_gain
    h_int
    function WaveSet(filename, hv)
        waveforms_raw = h5read(filename, "$hv/waveforms")
        v_gain = h5read(filename, "$hv/waveform_info/v_gain")
        h_int = h5read(filename, "$hv/waveform_info/h_int")
        new(waveforms_raw .* v_gain, v_gain, h_int)
    end
end

"""
$(SIGNATURES)
calculates charges of waveforms of a waveset 
# Arguments
- `waveforms => Matrix{Float64}`: waveforms
- `ped_min`: start of pedestal integration window
- `ped_max`: end of pedestal integration window
- `sig_min`: start of signal integration window
- `sig_min`: end of signal integration window
"""
function calculate_charges(waveforms, ped_min, ped_max, sig_min, sig_max)
    ped_sig_ratio = (ped_max - ped_min) / (sig_max - sig_min)
    pedestal = sum(waveforms[ped_min:ped_max, :], dims=1)
    charges = -(sum(waveforms[sig_min:sig_max, :], dims=1) - pedestal / ped_sig_ratio)
    vec(charges);
end

function bin_data(data, steprange::StepRangeLen)
    hist = append!(Histogram(steprange), data)
    edges = collect(hist.edges[1])
    x = edges .+ (edges[2] - edges[1]) / 2
    x = x[1:end-1]
    y = hist.weights
    x, y
end

struct ChargeDist
    x
    y
    """
    $(SIGNATURES)
    constructor
    # Arguments
    - `charges => Vector{Float64}`: charges from which charge distribution is 
                                    constructed
    - `steprange => StepRangeLen`: steprange for binning of the charges
    """
    function ChargeDist(charges, steprange::StepRangeLen)
        x, y = bin_data(charges, steprange)
        new(x, y)
    end
    
end
"""
    $(SIGNATURES)
    alternative constructor
    # Arguments
    - `charges => Vector{Float64}`: charges from which charge distribution is 
                                    constructed
    - `bins => Integer`: number of bins for binning of the charges
    - `lower`: lower edge for binning of the charges; default: -1
    - `upper`: upper edge for binning of the charges; default: 1
    """
ChargeDist(charges, bins::Integer; lower=-1, upper=1) = ChargeDist(charges, range(lower, upper; length=bins))


"""

abstract type FitResults end


struct PreFit <: FitResults
    μₚ
    σₚ
    Aₚ
    μₛ
    σₛ
    Aₛ
end

struct PMTRespFit <: FitResults
    μₚ
    σₚ
    μₛ
    σₛ
    nₚₑ
    A 
end

struct PMTRespUapFit <: FitResults
    μₚ
    σₚ
    μₛ
    σₛ
    nₚₑ
    A
    μᵤₐₚ
    σᵤₐₚ
    Aᵤₐₚ  
end


    $(SIGNATURES)
    normalized gaussian function
"""
norm_gauss(x, p) = exp.(.-(x .- p[1]).^2/(2 .* p[2] .^ 2)) ./ sqrt.(2π) ./ p[2]


"""
    $(SIGNATURES)
    gaussian function
"""
gauss(x, p) = p[3] .* exp.(.-(x .- p[1]).^2/(2 .* p[2] .^ 2)) ./ sqrt.(2π) ./ p[2]


"""
    $(SIGNATURES)
    PMT response function
"""
function pmtresp(x, p)
    poisson = Poisson(abs(p[5]))
    n_gaussians = 10
    pedestal = pdf(poisson, 0) .* norm_gauss(x, [p[1], p[2]])
    signal = zeros(length(x))
    for i in 1:n_gaussians
        signal .+= pdf(poisson, i) .* norm_gauss(x, [p[3] .* i, sqrt(i) .* p[4]])
    end
    p[6] .* (pedestal .+ signal)
end

function pmtresp(x, fit::PMTRespFit)
    pmtresp(x, [fit.μₚ, fit.σₚ, fit.μₛ, fit.σₛ, fit.nₚₑ, fit.A])
end


"""
    $(SIGNATURES)
    PMT response function with under amplified pulses modification
"""
function pmtresp_uap(x, p)
    poisson = Poisson(abs(p[5]))
    n_gaussians = 10
    pedestal = pdf(poisson, 0) .* norm_gauss(x, [p[1], p[2]])
    signal = zeros(length(x))
    for i in 1:n_gaussians
        signal .+= pdf(poisson, i) .* norm_gauss(x, [p[3] .* i, sqrt(i) .* p[4]])
    end
    p[6] .* (pedestal .+ signal) .+ gauss(x, [p[7], p[8], p[9]])
end

function pmtresp_uap(x, fit::PMTRespUapFit)
    pmtresp_uap(x, [fit.μₚ, fit.σₚ, fit.μₛ, fit.σₛ, fit.nₚₑ, fit.A, fit.μᵤₐₚ, fit.σᵤₐₚ, fit.Aᵤₐₚ])
end



function baseline_max(waveforms)
    floor(Int32, argmin(mean(waveforms, dims=2)).I[1] * 0.75)
end

function calculate_transit_times(waveforms, threshold)
    waveforms = waveforms .- mean(waveforms[1:baseline_max(waveforms), :], dims=1)
    n, m = size(waveforms)
    transit_times = zeros(m)
    for j in 1:m
        for i in 1:n
            value = waveforms[i,j]
            if value < threshold
                transit_times[j] = i - 1
                break
            end
        end
    end
    transit_times[transit_times .!= 0]
end