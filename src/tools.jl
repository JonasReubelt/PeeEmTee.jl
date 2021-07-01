using HDF5

function high_voltages(filename)
    """
    reads high voltages of a hdf5 file with wavesets

    Parameters
    ----------
    filename: string
    
    Returns
    -------
    high_voltages: Vector{String}
    """
    h5 = h5open(filename)
    high_voltages = keys(h5)
    close(h5)
    high_voltages
end

