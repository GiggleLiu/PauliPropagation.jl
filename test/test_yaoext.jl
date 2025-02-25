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

@testset "simulation" begin
    nq = 16
    pstr = PauliString(nq, :Z, 12)
    nl = 4
    topology = bricklayertopology(nq; periodic=false)

    circuit = hardwareefficientcircuit(nq, nl; topology=topology)
    nparams = countparameters(circuit)
    thetas = randn(nparams) * 0.5;

    yao_circ = yao_circuit(nq, circuit, thetas)
    psum_exact = propagate(circuit, pstr, thetas; min_abs_coeff=0)
    psum_yao = apply!(zero_state(nq), yao_circ)

    exact_expectation = overlapwithzero(psum_exact)
    yao_expectation = expect(put(nq, 12 => Z), psum_yao)
    @test isapprox(exact_expectation, yao_expectation, atol=1e-6)
end