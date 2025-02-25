"""
    yao_circuit(n::Int, circ::AbstractVector{Gate}, thetas::AbstractVector)

Convert a a circuit represented by a vector of gates and parameters to a Yao circuit representation that suited for exact simulation.

# Arguments
- `n::Int`: Number of qubits.
- `circ::AbstractVector{Gate}`: Vector of gates.
- `thetas::AbstractVector`: Vector of parameters.
"""
function yao_circuit(args...)
    error("You must `using Yao` first to use this feature.")
end

"""
    from_yao(circ::Chain)

Convert a Yao circuit to a vector of gates and parameters.

# Arguments
- `circ::Chain`: Yao circuit in the form of chain of basic gates.
"""
function from_yao(args...)
    error("You must `using Yao` first to use this feature.")
end