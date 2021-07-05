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
        new(waveforms_raw * v_gain, v_gain, h_int)
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
        hist = append!(Histogram(steprange), charges)
        edges = collect(hist.edges[1])
        x = edges .+ (edges[2] - edges[1]) / 2
        x = x[1:end-1]
        y = hist.weights
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

