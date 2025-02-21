using PauliPropagation

# number of qubits
nq = 16

# Pauli string with Z on qubit 32
pstr = PauliString(nq, :Z, 12)

# Number of layers
nl = 4

# Topology
topology = bricklayertopology(nq; periodic=false)

# Hardware-efficient circuit
circuit = hardwareefficientcircuit(nq, nl; topology=topology)

using Yao
function yao_circuit(n::Int, circ::Vector{Gate}, thetas::AbstractVector)
    @assert length(thetas) == countparameters(circuit)
    thetas = copy(thetas)
    c = chain(n)
    for g in circ
        if g isa PauliRotation
            push!(c, put(n, (g.qinds...,) => rot(kron(symbol_to_pauli.(g.symbols)...), popfirst!(thetas))))
        else
            error("Unsupported gate type: $(typeof(g))")
        end
    end
    return c
end
symbol_to_pauli(sym::Symbol) = sym == :X ? Yao.X : sym == :Y ? Yao.Y : Yao.Z

using Random
Random.seed!(42)
nparams = countparameters(circuit)
thetas = randn(nparams) * 0.5;


yao_circ = yao_circuit(nq, circuit, thetas)

@time psum_exact = propagate(circuit, pstr, thetas; min_abs_coeff=0)
@time psum_yao = apply!(zero_state(nq), yao_circ)

exact_expectation = overlapwithzero(psum_exact)
yao_expectation = expect(put(nq, 12 => Z), psum_yao)