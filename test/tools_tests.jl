using PeeEmTee
using Test
using StatsBase

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

@testset "calculate_transit_times" begin
    a = [[0., -1., 1., -53., 0., 1.] [-5., 3., -14., -253., 143., 3.] [100., 102., 85., -9., -1., -15.]]
    @test calculate_transit_times(a, -10.) == [3, 3, 2]
end

@testset "bin_data()" begin
    a = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    @test bin_data(a, 1.:1.:10.) == ([1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5],
                                     [1, 1, 1, 1, 1, 1, 1, 1, 1])
end

@testset "ChargeDist" begin
    a = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    chargedist = ChargeDist(a, 1.:1.:10.)
    @test chargedist.x == [1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5]
    @test chargedist.y == [1, 1, 1, 1, 1, 1, 1, 1, 1]
    b = [-.8, -0.3, 0.1, 0.5]
    chargedist = ChargeDist(b, 4)
    @test chargedist.x == [-0.75, -0.25, 0.25, 0.75]
    @test chargedist.y == [1, 1, 1, 1]
end

@testset "norm_gauss()" begin
    @test norm_gauss(0, [0, 1]) == 0.3989422804014327
end

@testset "gauss()" begin
    @test gauss(0, [0, 1, 10]) == 3.9894228040143274
end

@testset "baseline_max" begin
    a = [[0., 0., -1., 0., 0.] [0., 0., -2., 0., -1.] [0., 1., -3., -2., -1.]]
    @test baseline_max(a) == 2
end

@testset "subtract_baseline()" begin
    a = [[-1., 0., -2., 0., 0.] [0., 0., -2., 0., -1.] [0., 1., -3., -2., -1.]]
    a_res1 = [[-0.5, 0.5, -1.5, 0.5, 0.5] [0., 0., -2., 0., -1.] [-0.5, 0.5, -3.5, -2.5, -1.5]]
    a_res2 = [[0, 1., -1., 1., 1.] [0., 0., -2., 0., -1.] [0., 1., -3., -2., -1.]]
    @test subtract_baseline(a) == a_res1
    @test subtract_baseline(a; baseline_range=(1, 1)) == a_res2
end

@testset "pmtresp()" begin
    xs = [0, 0.3]
    params = [0, 0.05, 0.3, 0.2, 0.1, 100.]
    @test pmtresp(xs, params) == [727.8832805768451,
                                  18.4164040244264]
    prffit = PMTRespFit(params...)
    @test pmtresp(xs, prffit) == [727.8832805768451,
                                  18.4164040244264]
end

@testset "pmtresp_uap()" begin
    xs = [0.1, 0.2]
    params = [0, 0.05, 0.3, 0.2, 0.1, 100., 0.1, 0.02, 10.]
    @test pmtresp_uap(xs, params) == [308.25941556084206,
                                      16.40805075212658]
    prffit = PMTRespUapFit(params...)
    @test pmtresp_uap(xs, prffit) == [308.25941556084206,
                                      16.40805075212658]
end

@testset "simulate_charges()" begin
    secs = (7, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5)
    nₚₑ = 0.1
    @test isapprox(mean(simulate_charges(nₚₑ, 100000; secs=secs)),
                   prod(secs) * nₚₑ;
                   rtol=0.05) 
end