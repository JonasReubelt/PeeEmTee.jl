using PeeEmTee
using Test
using DelimitedFiles

const chargedist_file = joinpath(@__DIR__, "data", "chargedist_example.txt")

@testset "make_qfunc()" begin
    data = readdlm(chargedist_file)
    params = [-0.0009971235435351372,
              0.01266241070136633,
              0.1999583888727496,
              0.0996819367396509,
              0.09768461062509437,
              499.5820279060795]
    qfunc = make_qfunc(pmtresp, data[:, 1], data[:, 2])
    @test qfunc(params) == 66.02414426599655 
end

@testset "pre_fit()" begin
    data = readdlm(chargedist_file)
    chargedist = ChargeDist(data[:, 1], data[:, 2])
    prefit = PreFit(-0.0009439575776608711,
                    0.012781642902921997,
                    455.31645300021285,
                    0.20214663928517476,
                    0.10553835465661239,
                    46.63963250021425)
    @test pre_fit(chargedist) == prefit
end