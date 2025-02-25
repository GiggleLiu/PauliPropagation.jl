module YaoExt

using Yao, PauliPropagation

function PauliPropagation.yao_circuit(n::Int, circ::AbstractVector{Gate}, thetas::AbstractVector)
    @assert length(thetas) == countparameters(circ)
    thetas = copy(thetas)
    c = chain(n)
    for g in circ
        if g isa PauliRotation
            push!(c, put(n, (g.qinds...,) => rot(kron(symbol_to_yao.(g.symbols)...), popfirst!(thetas))))
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
