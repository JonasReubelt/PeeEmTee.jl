using PeeEmTee
using Test
using DelimitedFiles

const chargedist = joinpath(@__DIR__, "data", "chargedist_example.txt")

@testset "make_qfunc()" begin
    data = readdlm(chargedist)
    params = [-0.0009971235435351372,
              0.01266241070136633,
              0.1999583888727496,
              0.0996819367396509,
              0.09768461062509437,
              499.5820279060795]
    qfunc = make_qfunc(pmtresp, data[:, 1], data[:, 2])
    @test qfunc(params) == 66.02414426599655 
end