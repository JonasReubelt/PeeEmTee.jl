using PeeEmTee
using Test

const pmt_data = joinpath(@__DIR__, "data", "pmt_data.h5")


@testset "WaveSet" begin
    a = [[1, 2, 3] [4, 5, 6] [7, 8, 9]]
    waveset = WaveSet(pmt_data, 1000)
    @test waveset.waveforms == a * waveset.v_gain
    @test waveset.h_int == 1.0
    waveset = WaveSet(pmt_data, 1100)
    @test waveset.waveforms == a * waveset.v_gain
    @test waveset.h_int == 1.0
    waveset = WaveSet(pmt_data, 1200)
    @test waveset.waveforms == a * waveset.v_gain
    @test waveset.h_int == 1.0
end