using Yao, PauliPropagation, Test
using Random

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


# noisy simulation
@testset "pauli noisy simulation" begin
    Random.seed!(1234)
    nq = 3
    nl = 2
    W = Inf
    min_abs_coeff = 0.0
    opind = rand(1:nq)
    pstr = PauliString(nq, :Z, opind)

    topo = bricklayertopology(nq; periodic=false)
    circ = hardwareefficientcircuit(nq, nl; topology=topo)

    m = countparameters(circ)

    depolarizing_circ = deepcopy(circ)
    pauli_circ = deepcopy(circ)

    where_ind = rand(1:m)
    q_ind = opind
    noise_p = 0.1 # rand() * 0.2
    insert!(depolarizing_circ, where_ind, DepolarizingNoise(q_ind))
    insert!(pauli_circ, where_ind, PauliZNoise(q_ind))
    insert!(pauli_circ, where_ind, PauliYNoise(q_ind))
    insert!(pauli_circ, where_ind, PauliXNoise(q_ind))

    thetas1 = rand(m)
    thetas2 = deepcopy(thetas1)
    thetas3 = deepcopy(thetas1)
    insert!(thetas1, where_ind, noise_p)
    insert!(thetas2, where_ind, noise_p)
    insert!(thetas2, where_ind, noise_p)
    insert!(thetas2, where_ind, noise_p)
    pauli_p = 0.5 - sqrt(0.25 - noise_p/4)
    insert!(thetas3, where_ind, pauli_p)
    insert!(thetas3, where_ind, pauli_p)
    insert!(thetas3, where_ind, pauli_p)


    dnum1 = propagate(depolarizing_circ, pstr, thetas1; max_weight=W, min_abs_coeff=min_abs_coeff)
    dnum2 = propagate(pauli_circ, pstr, thetas2; max_weight=W, min_abs_coeff=min_abs_coeff)

    exp1 = overlapwithzero(dnum1)
    exp2 = overlapwithzero(dnum2)
    @test exp1 ≈ exp2

    # convert to Yao blocks
    yao_depolarizing_circ = yao_circuit(nq, depolarizing_circ, thetas1)
    yao_pauli_circ = yao_circuit(nq, pauli_circ, thetas3)
    exp_yao1 = expect(put(nq, opind => Z), zero_state(nq) |> density_matrix |> yao_depolarizing_circ)
    exp_yao2 = expect(put(nq, opind => Z), zero_state(nq) |> density_matrix |> yao_pauli_circ)

    @test exp1 ≈ exp_yao1
    @test exp2 ≈ exp_yao2

    # convert back
    n, circ, thetas = from_yao(yao_depolarizing_circ; frozen_rots=false)
    @test circ == depolarizing_circ
    @test thetas == thetas1
    n, circ, thetas = from_yao(yao_pauli_circ; frozen_rots=false)
    @test circ == pauli_circ
    @test thetas == thetas3
end

@testset "dephasing noisy simulation" begin
    nq = 16
    nl = 4
    W = Inf
    min_abs_coeff = 0.0
    opind = rand(1:nq)
    pstr = PauliString(nq, :Z, opind)

    topo = bricklayertopology(nq; periodic=false)
    circ = hardwareefficientcircuit(nq, nl; topology=topo)

    m = countparameters(circ)

    dephasing_circ = deepcopy(circ)
    pauli_circ = deepcopy(circ)

    where_ind = rand(1:m)
    q_ind = opind
    noise_p = randn() * 0.2
    insert!(dephasing_circ, where_ind, DephasingNoise(q_ind))
    insert!(pauli_circ, where_ind, PauliYNoise(q_ind))
    insert!(pauli_circ, where_ind, PauliXNoise(q_ind))

    thetas1 = randn(m)
    thetas2 = deepcopy(thetas1)
    insert!(thetas1, where_ind, noise_p)
    insert!(thetas2, where_ind, noise_p)
    insert!(thetas2, where_ind, noise_p)

    dnum1 = propagate(dephasing_circ, pstr, thetas1; max_weight=W, min_abs_coeff=min_abs_coeff)

    dnum2 = propagate(pauli_circ, pstr, thetas2; max_weight=W, min_abs_coeff=min_abs_coeff)

    @test overlapwithzero(dnum1) ≈ overlapwithzero(dnum2)
end