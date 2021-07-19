using PeeEmTee
using Test

const pmt_data = joinpath(@__DIR__, "data", "pmt_data.h5")

@testset "high_voltages()" begin
    @test high_voltages(pmt_data) == ["1000", "1100", "1200"]
end

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

@testset "calculate_charges" begin
    a = [[0, 1, -45, -53, 0, -1] [-5, 3, -145, -253, 3, -5] [0, 0, -44, -12, -1, 1]]
    @test calculate_charges(a, 1, 2, 3, 4) == [99, 396, 56]
end