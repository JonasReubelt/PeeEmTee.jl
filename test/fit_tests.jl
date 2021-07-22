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

@testset "pmtresp_fit()" begin
    data = readdlm(chargedist_file)
    chargedist = ChargeDist(data[:, 1], data[:, 2])
    prffit = PMTRespFit(-0.0009971235435351222,
                        0.0126624107013663,
                        0.1999583888727531,
                        0.09968193673964637,
                        0.0976846106250966,
                        499.58202790610073)
    prefit = pre_fit(chargedist)
    @test pmtresp_fit(chargedist, prefit) == prffit
    prfuapfit = PMTRespUapFit(-0.0010866013489182497,
                              0.012573335066167545,
                              0.20255116463024672,
                              0.09761207500313138,
                              0.09719044076647844,
                              495.05266162357447,
                              0.012985248350779485,
                              0.021308630425982703,
                              4.592234142351302)
    @test pmtresp_fit(chargedist, prefit, mod=:uap) == prfuapfit
end