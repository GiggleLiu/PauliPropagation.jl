using Yao, PauliPropagation, Test

@testset "Test YaoExt" begin
    nq = 3
    nl = 2
    circuit = tfitrottercircuit(nq, nl)
    yaocirc = yao_circuit(nq, circuit, randn(countparameters(circuit)))
    n, circ, thetas = from_yao(yaocirc; frozen_rots=false)
    @test n == nq
    @test circ == circuit
    @test thetas == thetas
end