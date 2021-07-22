"""
$(SIGNATURES)
reads the high voltages of the wavesets saved in the hdf5 file
# Arguments
- `filename`: filename of h5 file
"""
function high_voltages(filename)
    h5 = h5open(filename, "r")
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
end

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
        ChargeDist(x, y)
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
ChargeDist(charges, bins::Integer; lower=-1, upper=1) = ChargeDist(charges, range(lower, upper; length=bins+1))




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

"""
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


"""
    $(SIGNATURES)
    calculates maximum for baseline calculation of a set of waveforms
    # Arguments
    - `waveforms => Matrix{Float64}`: set of waveforms
"""
function baseline_max(waveforms::Matrix{Float64})
    floor(Int32, argmin(mean(waveforms, dims=2)).I[1] * 0.75)
end

"""
    $(SIGNATURES)
    subtracts baseline from a set of waveforms, i.e. sets baseline to 0
    # Arguments
    - `waveforms => Matrix{Float64}`: set of waveforms
    - `baseline_range => Tuple`: range in which mean of baseline is calculated
                                 (default: nothing => baseline_max is used to calculate range)
"""
function subtract_baseline(waveforms::Matrix{Float64}; baseline_range=nothing)
    if baseline_range==nothing
        baseline_range = (1, baseline_max(waveforms))
    end
    waveforms .- mean(waveforms[baseline_range[1]:baseline_range[2], :], dims=1)
end

"""
    $(SIGNATURES)
    calculates transit times of a set of waveforms
    # Arguments
    - `waveforms => Matrix{Float64}`: set of waveforms
    - `threshold`: first pass of threshold => transit time
"""
function calculate_transit_times(waveforms::Matrix{Float64}, threshold)
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

"""
    $(SIGNATURES)
    shifts waveforms by its minimum and calculates mean waveform
    # Arguments
    - `waveforms => Matrix{Float64}`: set of waveforms 
"""
function mean_waveform(waveforms::Matrix{Float64})
    a, b = size(waveforms)
    shifted_waveforms = Matrix{Float64}(undef, a, b)
    for j in 1:b
        minval, minidx = findmin(waveforms[:, j])
        for i in 1:a
            shift = floor(Int32, a/2)
            shifted_i = i - minidx + shift
            if shifted_i > a
                shifted_waveforms[shifted_i - a, j] = waveforms[i, j]
            elseif shifted_i < 1
                shifted_waveforms[shifted_i + a, j] = waveforms[i, j]
            else
                shifted_waveforms[shifted_i, j] = waveforms[i, j]
            end
        end
    end
    vec(mean(shifted_waveforms, dims=2))
end

"""
    $(SIGNATURES)
    simulations PMT charge spectrum with Poisson distributed
    secondary emission coefficients at each dynode
    # Arguments
    - `nₚₑ`: mean number of photoelectrons
    - `N`: number of simulated charges
    - `σₚ`: Gaussian sigma of pedestal (default: 1e5)
    - `secs`: secondary emission coefficients
              (default: (7, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5))
"""
function simulate_charges(nₚₑ, N;
                          σₚ=1e5, 
                          secs=(7, 4.5, 4.5, 4.5, 4.5,
                                4.5, 4.5, 4.5, 4.5, 4.5))
    charges = Float64[]
    for i in 1:N
        pes = rand(Poisson(nₚₑ))
        gaussian = Normal(0., σₚ)
        if pes == 0
            push!(charges, rand(gaussian))
        else
            charge = 0
            for i in 1:pes
                e = rand(Poisson(secs[1]))
                for i in 2:10
                    e = rand(Poisson(e * secs[i]))
                end
                charge += e
            end
            push!(charges, charge + rand(gaussian))
        end
    end
    charges
end

