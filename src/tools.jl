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


"""
    $(SIGNATURES)
    normalized gaussian function
"""
norm_gauss(x, p) = exp.(.-(x .- p[1]).^2/(2 .*p[2] .^ 2)) ./ sqrt.(2π) ./ p[2]


"""
    $(SIGNATURES)
    gaussian function
"""
gauss(x, p) = p[3] .* exp.(.-(x .- p[1]).^2/(2 .*p[2] .^ 2)) ./ sqrt.(2π) ./ p[2]


"""
    $(SIGNATURES)
    PMT response function
"""
function pmtresp(x, p)
    poisson = Poisson(abs(p[5]))
    n_gaussians = 10
    pedestal = pdf(poisson, 0) .* gauss.(x, p[1], p[2])
    signal = zeros(length(x))
    for i in 1:n_gaussians
        signal .+= pdf(poisson, i) .* gauss.(x, p[3] * i, sqrt(i) * p[4])
    end
    p[6] .* (pedestal .+ signal)
end

"""
    $(SIGNATURES)
    quality function constructor for least square optimization
"""
function make_qfunc(model, x, y)
    function qfunc(p)
        mask = (y .!= 0)
        model_res = model(x, p)
        sum(((model_res[mask] .- y[mask]).^2) ./ y[mask])
    end
    qfunc
end


struct PreFitResults
    μₚ
    σₚ
    Aₚ
    μₛ
    σₛ
    Aₛ
end


"""
    $(SIGNATURES)
    performs pre fit
"""
function pre_fit(chargedist::ChargeDist)
    x = chargedist.x
    y = chargedist.y
    qfunc = make_qfunc(gauss, x, y)
    mxval_ped, mxidx_ped = findmax(y)
    p0 = [x[mxidx_ped], 0.01, mxval_ped]
    fit_ped = optimize(qfunc, p0, NewtonTrustRegion())
    popt_ped = Optim.minimizer(fit_ped)
    
    mask = (x .> popt_ped[1] + 5 * popt_ped[2])
    x_spe = x[mask]
    y_spe = y[mask]
    qfunc = make_qfunc(gauss, x_spe, y_spe)
    mxval_spe, mxidx_spe = findmax(y_spe)
    p0 = [x_spe[mxidx_spe], 0.1, mxval_spe]
    fit_spe = optimize(qfunc, p0, NewtonTrustRegion())
    popt_spe = Optim.minimizer(fit_spe)
    PreFitResults(popt_ped..., popt_spe...)
end
