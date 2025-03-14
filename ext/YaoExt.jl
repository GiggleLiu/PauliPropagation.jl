module YaoExt

using Yao, PauliPropagation

function PauliPropagation.yao_circuit(n::Int, circ::AbstractVector{<:Gate}, thetas::AbstractVector)
    @assert length(thetas) == countparameters(circ)
    thetas = copy(thetas)
    c = chain(n)
    for g in circ
        if g isa PauliRotation
            push!(c, put(n, (g.qinds...,) => rot(kron(symbol_to_yao.(g.symbols)...), popfirst!(thetas))))
        elseif g isa DepolarizingNoise
            @assert length(g.qind) == 1 "Depolarizing noise should be applied to a single qubit"
            push!(c, put(n, (g.qind...,) => single_qubit_depolarizing_channel(popfirst!(thetas))))
        elseif g isa PauliXNoise
            push!(c, put(n, (g.qind...,) => pauli_error_channel(; px=popfirst!(thetas), py=0.0, pz=0.0)))
        elseif g isa PauliYNoise
            push!(c, put(n, (g.qind...,) => pauli_error_channel(; px=0.0, py=popfirst!(thetas), pz=0.0)))
        elseif g isa PauliZNoise
            push!(c, put(n, (g.qind...,) => pauli_error_channel(; px=0.0, py=0.0, pz=popfirst!(thetas))))
        elseif g isa DephasingNoise
            error("Dephasing noise not implemented")
        else
            error("Unsupported gate type: $(typeof(g))")
        end
    end
    return c
end

function PauliPropagation.from_yao(circ::ChainBlock; frozen_rots::Bool=false)
    circ = Yao.Optimise.to_basictypes(circ)
    n = nqubits(circ)
    gates = Gate[]
    thetas = Float64[]
    for g in circ
        if g isa PutBlock && g.content isa RotationGate
            if frozen_rots
                push!(gates, PauliRotation(yao_to_symbols(g.content.block), collect(g.locs), g.content.theta))
            else
                push!(gates, PauliRotation(yao_to_symbols(g.content.block), collect(g.locs)))
                push!(thetas, g.content.theta)
            end
        elseif g isa PutBlock && g.content isa UnitaryChannel
            for (prob, operator) in zip(g.content.probs, g.content.operators)
                if operator isa Yao.XGate
                    push!(gates, PauliYNoise(collect(g.locs), prob))
                    push!(gates, PauliZNoise(collect(g.locs), prob))
                elseif operator isa YGate
                    push!(gates, PauliZNoise(collect(g.locs), prob))
                    push!(gates, PauliXNoise(collect(g.locs), prob))
                elseif operator isa ZGate
                    push!(gates, PauliXNoise(collect(g.locs), prob))
                    push!(gates, PauliYNoise(collect(g.locs), prob))
                elseif !(operator isa Yao.I2Gate)
                    error("Unsupported error type: $(typeof(operator))")
                end
            end
        else
            error("Unsupported gate type: $(typeof(g))")
        end
    end
    return n, gates, thetas
end

const symbol_yao_map = Dict(
    :I => ConstGate.I2,
    :X => ConstGate.X,
    :Y => ConstGate.Y,
    :Z => ConstGate.Z,
    :H => ConstGate.H,
    :S => ConstGate.S,
    :T => ConstGate.T,
)
const yao_symbol_map = Dict(values(symbol_yao_map) .=> keys(symbol_yao_map))

function symbol_to_yao(sym::Symbol)
    return symbol_yao_map[sym]
end

function yao_to_symbols(yaoblock::KronBlock)
    return [yao_to_symbol(b) for b in yaoblock.blocks]
end
yao_to_symbols(yaoblock) = [yao_to_symbol(yaoblock)]

function yao_to_symbol(yaoblock)
    return yao_symbol_map[yaoblock]
end

end
