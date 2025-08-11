using PauliPropagation

# number of qubits
nq = 64

# Pauli string with Z on qubit 32
pstr = PauliString(nq, :Z, 32)
# Example: PauliString(4, [:X, :Y, :Z], [2, 3, 4], 2.0)
# input_psum = PauliSum(nq)
# add!(input_psum, :Z, 32)
# add!(input_psum, [:X, :X], [1, 64])
# add!(input_psum, [:X, :Y, :Z], [1, 4, 7], 1.3)

# Number of layers
nl = 4

# Topology
topology = bricklayertopology(nq; periodic=false)

# Hardware-efficient circuit
circuit = hardwareefficientcircuit(nq, nl; topology=topology)

# Number of parameters
nparams = countparameters(circuit)

# rectopology = rectangletopology(8, 8; periodic=false) 

# hardcoded_custom_topology = [(1, 2), (1, 15), (2, 3), (3, 4), (4, 5), (5, 6), (5, 16), (6, 7)]

# function centralspintopology(nqubits)
#     return [(1, ii) for ii in 2:nqubits]
# end

# centralspin_topology = centralspintopology(nq)

# circuit = tfitrottercircuit(nq, nl; topology=topology, start_with_ZZ=true);
# circuit = tiltedtfitrottercircuit(nq, nl; topology=topology)
# circuit = heisenbergtrottercircuit(nq, nl; topology=topology)

# function myheisenbergtrottercircuit(nqubits::Integer, nlayers::Integer; topology=nothing)
#     circuit::Vector{Gate} = []

#     if isnothing(topology)
#         topology = bricklayertopology(nqubits)
#     end

#     for _ in 1:nlayers
#         rxxlayer!(circuit, topology)
#         ryylayer!(circuit, topology)
#         rzzlayer!(circuit, topology)
#     end

#     return circuit
# end

# circuit = myheisenbergtrottercircuit(nq, nl; topology=topology);
# nparams = countparameters(circuit)

# Jxx = 0.5
# Jyy = 0.25
# Jzz = 0.25
# dt = 0.2

# thetaxx_indices = getparameterindices(circuit, PauliRotation, [:X, :X])
# thetayy_indices = getparameterindices(circuit, PauliRotation, [:Y, :Y])
# thetazz_indices = getparameterindices(circuit, PauliRotation, [:Z, :Z])

# thetas = zeros(nparams)
# thetas[thetaxx_indices] .= 2 * Jxx * dt
# thetas[thetayy_indices] .= 2 * Jyy * dt
# thetas[thetazz_indices] .= 2 * Jzz * dt;

using Random
Random.seed!(42)
thetas = randn(nparams) * 0.5;

@time psum = propagate(circuit, pstr, thetas)

@time propagate(circuit, pstr, thetas);
overlapwithzero(psum)

@time psum_exact = propagate(circuit, pstr, thetas; min_abs_coeff=0)
print("The error with our default truncation of `1e-10` is ", abs(overlapwithzero(psum_exact) - overlapwithzero(psum)))

min_abs_coeff = 1e-3
@time psum_coeff = propagate(circuit, pstr, thetas; min_abs_coeff = min_abs_coeff)

print("The error with `min_abs_coeff = $min_abs_coeff` is ", abs(overlapwithzero(psum_exact) - overlapwithzero(psum_coeff)))

# input_psum = PauliSum(nq)
# add!(input_psum, :Z, 32)
# add!(input_psum, [:X, :X], [1, 64])
# add!(input_psum, [:X, :Y, :Z], [1, 4, 7], 1.3)

# @time propagate!(circuit, input_psum, thetas; min_abs_coeff=1e-3)
# @time propagate!(circuit, input_psum, thetas; min_abs_coeff=1e-3) # and again
# input_psum

max_weight = 5
@time psum_weight = propagate(circuit, pstr, thetas; max_weight = max_weight)

print("The error with `max_weight = $max_weight` is ", abs(overlapwithzero(psum_exact) - overlapwithzero(psum_weight)))

trunc_coeffs = 10.0 .^ (-1:-1:-10)
# One more exact computation. Record the expectation value and time.
res = @timed propagate(circuit, pstr, thetas; min_abs_coeff=0)
exact_expectation = overlapwithzero(res.value)
exact_time = res.time

# Sweep over the truncation values
expectations = zeros(length(trunc_coeffs))
times = zeros(length(trunc_coeffs))
for (ii, min_abs_coeff) in enumerate(trunc_coeffs)
    res = @timed propagate(circuit, pstr, thetas; min_abs_coeff=min_abs_coeff)
    expectations[ii] = overlapwithzero(res.value)
    times[ii] = res.time
end

# calculate the absolute error to the exact expectation value
abs_errors = abs.(expectations .- exact_expectation);


using CairoMakie

# Create a figure with 4 subplots arranged in a 2x2 grid
fig = Figure(size=(1000, 800))

# Plot 1: Error vs Truncation Threshold

ax1 = Axis(fig[1,1], 
    xscale=log10, 
    yscale=log10,
    xlabel="Truncation Threshold",
    ylabel="Error",
    xticks=trunc_coeffs,
    yticks=10.0 .^ (-10:1:0))
scatter!(ax1, trunc_coeffs, abs_errors)
lines!(ax1, trunc_coeffs, abs_errors)

# Plot 2: Runtime vs Truncation Threshold  
ax2 = Axis(fig[1,2],
    xscale=log10,
    yscale=log10, 
    xlabel="Truncation Threshold",
    ylabel="Runtime [s]",
    xticks=trunc_coeffs,
    yticks=10.0 .^ (-4:1:0))
scatter!(ax2, trunc_coeffs, times)
lines!(ax2, trunc_coeffs, times, label="With coefficient truncation")
hlines!(ax2, [exact_time], color=:grey, linestyle=:dash, linewidth=2, label="Exact computation")
axislegend(ax2)

# Plot 3: Runtime vs Error
ax3 = Axis(fig[2,1],
    xscale=log10,
    yscale=log10,
    xlabel="Error",
    ylabel="Runtime [s]",
    xticks=10.0 .^ (-10:2:0),
    yticks=10.0 .^ (-4:1:0))
scatter!(ax3, abs_errors, times, label="With coefficient truncation")
lines!(ax3, abs_errors, times)
hlines!(ax3, [exact_time], color=:grey, linestyle=:dash, linewidth=2, label="Exact computation")
axislegend(ax3)

# Plot 4: Expectation Value vs Truncation Threshold
ax4 = Axis(fig[2,2],
    xscale=log10,
    xlabel="Truncation Threshold",
    ylabel="Expectation Value",
    xticks=trunc_coeffs)
scatter!(ax4, trunc_coeffs, expectations)
lines!(ax4, trunc_coeffs, expectations)

# Adjust spacing between subplots
fig[1:2,1:2] = GridLayout()

# Display the figure
fig