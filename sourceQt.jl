using LinearAlgebra

function density_matrix(state::state_vector)
    return density_matrix(state.vec*state.vec', state.base)
end

function pauli_measure(op::pauli_operator, qubit::Int64, nqubit::Int64)
    A = op.evec[1]*op.evec[1]'
    B = op.evec[2]*op.evec[2]'
    X = 1
    Y = 1
    for i in 1:nqubit
        if i == qubit
            X = kron(X, A)
            Y = kron(Y, B)
        else
            X = kron(X, Matrix(I, 2, 2))
            Y = kron(Y, Matrix(I, 2, 2))
        end
    end
    return pauli_measure([X, Y], op.base, op.eval)
end

#=function two_qubit_simulation_tree(state, measurements, channel)
    ρ_in = density_matrix(state)
    root = create_tree(measurements[1:2], ρ_in)
    leaves = find_leaves(root)
    for leaf in leaves
        ρ_out = density_matrix(channel.mat*leaf.value.mat*adjoint(channel.mat), ρ_in.base)
        attach_node!(root, create_tree(measurements[3:4], ρ_out, start = leaf.idx, parent = leaf.parent, prob=leaf.value.prob, res=leaf.value.result))
    end
    return root
end
=#

#=to re-tink
function measures(root::TreeNode, A)
    if root.left === nothing && root.right === nothing
        return
    end
    if root.value.result == 0
        return
    end
    rand() < root.left.value.prob ? (push!(A, root.left.value.result), measures(root.left, A)) : (push!(A, root.right.value.result), measures(root.right, A))
end
=#

function measures(root::TreeNode, A)
    if root.left === nothing && root.right === nothing
        return
    end
    if root.left.value.result == 0
        measures(root.left, A)
    end
    if root.right.value.result == 0
        measures(root.right, A)
    end
    rand() < root.left.value.prob ? (push!(A, root.left.value.result), measures(root.left, A)) : (push!(A, root.right.value.result), measures(root.right, A))
end

function simulation(Qt::QT_Dynamical_model, n)
    ρ_in = density_matrix(Qt.state)
    root = create_tree(Qt.operators, ρ_in)
    A = Vector{Vector{Int64}}(undef, n)
    B = Vector{Vector{Int64}}(undef, n)
    for i in 1:n
        G = []
        measures(root, G)
        A[i] = [G[1], G[3]]
        B[i] = [G[2], G[4]]
        #push!(A, [G[1], G[3]])
        #push!(B,[G[2], G[4]])
    end
    return A, B
end

σx = pauli_operator([0 1; 1 0], 1, "Z", [[1; 1]./sqrt(2), [-1; 1]./sqrt(2)], [1; -1])
σy = pauli_operator([0 -im; im 0], 2, "Z", [[-im; 1]./sqrt(2), [1;-im]./sqrt(2)], [1; -1])
σz = pauli_operator([1 0; 0 -1], 3, "Z", [[1; 0], [0; 1]], [1; -1])
Id = pauli_operator(Matrix(I, 2, 2), 0, "Z", [[1; 0], [0; 1]], [1]);
qubit = [1, 2]
nqubit = length(qubit)
Mx = [pauli_measure(σx, i, nqubit) for i in qubit]
My = [pauli_measure(σy, i, nqubit) for i in qubit]
Mz = [pauli_measure(σz, i, nqubit) for i in qubit]
MId = [pauli_measure(Id, i, nqubit) for i in qubit];


# Define quantum gates for 2-qubit systems
function cnotAB()
    return [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0]
end

function hadamard()
    return (1 / sqrt(2)) * [1 1; 1 -1]
end

function pauliX()
    return σx.mat
end

function swap()
    return [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1]
end

general_matrix(a,x,y,z) = Matrix{ComplexF64}(I, 2,2)*a + σx.mat*x + σy.mat*y + σz.mat*z

z_rotation(θ) = [cos(θ) -sin(θ) 0;
                  sin(θ) cos(θ) 0;
                  0 0 1]

x_rotation(ϕ) = [1 0 0;
                0 cos(ϕ) -sin(ϕ);
                0 sin(ϕ) cos(ϕ)]

rot(θ,ϕ) = z_rotation(θ)*x_rotation(ϕ)




gates = Dict(
    "CNOT_AB" => cnotAB(),
    "Hadamard" => hadamard(),
    "PauliX" => pauliX(),
    "Swap" => swap()
)

pauli_proj_meas = Dict(
    "Mx" => Dict("A" => Mx[1],
                 "B" => Mx[2]),
    "My" => Dict("A" => My[1],
                 "B" => My[2]),
    "Mz" => Dict("A" => Mz[1],
                 "B" => Mz[2]),
    "MId" => Dict("A" => MId[1],
                  "B" => MId[2])
)
