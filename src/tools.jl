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

