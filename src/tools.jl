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
    pedestal = sum(data[ped_min:ped_max, :], dims=1)
    charges = -(sum(data[sig_min:sig_max, :], dims=1) - pedestal / ped_sig_ratio)
    vec(charges);
end

