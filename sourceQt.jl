using LinearAlgebra

function density_matrix(state::state_vector)
    return density_matrix(state.vec*state.vec', state.base)
end

function noisy_density_matrix(state::state_vector, noise::Float64)
    return density_matrix(((1-noise)*(state.vec*state.vec') + (noise/4)*Matrix(I,4,4)), state.base)
end

function measure(op::operators, qubit::Int64, nqubit::Int64)
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
    return measure([X, Y], op.base, op.eval)
end

function operator(matrix, base)
    evals, evecs = eigen(matrix)
    eva = round.(evals)
    eve = [evecs[:,i] for i in 1:size(evecs)[2]]
    return operator(matrix, base, eve, eva)
end

function pauli_operator(matrix,base)
    evals, evecs = eigen(matrix)
    eve = [evecs[:,i] for i in 1:size(evecs)[2]]
    return pauli_operator(matrix, base, eve, evals)
end

function measures(root::TreeNode)
    if root.children[1] === nothing && root.children[2] === nothing
        return []
    end
    if root.children[1] === nothing && root.children[2] !== nothing
        if root.children[2].value.result == 0
            return [measures(root.children[2])...]
        else
            return [root.children[2].value.result; measures(root.children[2])...]
        end
        return
    end
    if root.children[2] === nothing && root.children[1] !== nothing
        if root.children[1].value.result == 0
            
            return [measures(root.children[1])...]
        else
            return [root.children[1].value.result; measures(root.children[1])...]
        end
    end
    if root.children[1].value.result == 0
        
        return [measures(root.children[1])...]
    end
    if root.children[2].value.result == 0
        return [measures(root.children[2])...]
    end
    rand() < root.children[1].value.prob ? (return[root.children[1].value.result; measures(root.children[1])...]) : (return[root.children[2].value.result; measures(root.children[2])...])
end
function noisy_simulation(Qt::QT_Dynamical_model, n, noise)
    ρ_in = noisy_density_matrix(Qt.state, noise)
    root = create_tree(Qt.operators, ρ_in, dummyBinary_tree())
    trimming!(root)
    A = Vector{Vector{Float64}}(undef, n)
    B = Vector{Vector{Float64}}(undef, n)
    for i in 1:n
        G = measures(root)
        A[i] = [G[1], G[3]]
        B[i] = [G[2], G[4]]
    end
    return A, B
end

function simulation(Qt::QT_Dynamical_model, n)
    ρ_in = density_matrix(Qt.state)
    root = create_tree(Qt.operators, ρ_in, dummyBinary_tree())
    trimming!(root)
    A = Vector{Vector{Float64}}(undef, n)
    B = Vector{Vector{Float64}}(undef, n)
    for i in 1:n
        G = measures(root)
        A[i] = [G[1], G[3]]
        B[i] = [G[2], G[4]]
    end
    return A, B
end

function rot(rotated::operators,rotator::pauli_operator, θ)
    A = (cos(θ/2)*Matrix(I, 2, 2) - im*sin(θ/2)*rotator.mat)*rotated.mat*(cos(θ/2)*Matrix(I, 2, 2) + im*sin(θ/2)*rotator.mat)
    B = real.(A)
    C = imag.(A)
    B[abs.(B) .< 1e-10] .= 0
    C[abs.(C) .< 1e-10] .= 0
    return B + im*C
end

function spherical_point(θ,ϕ)
    A = sin(θ)*cos(ϕ)*σx.mat + sin(θ)*sin(ϕ)*σy.mat + cos(θ)*σz.mat
    B = real.(A)
    C = imag.(A)
    B[abs.(B) .< 1e-10] .= 0
    C[abs.(C) .< 1e-10] .= 0
    return B + im*C
end

function projector(θ, ϕ, eig, qubit)
    if qubit == 1
        return kron((Matrix(I,2,2) +eig.*spherical_point(θ,ϕ))./2,Matrix(I,2,2))
    end
    if qubit == 2
        return kron(Matrix(I,2,2),(Matrix(I,2,2) +eig.*spherical_point(θ,ϕ))./2)
    end
    if qubit == 0
        return (Matrix(I,2,2) +eig.*spherical_point(θ,ϕ))./2
    end
end

function qtfullprob(pa,pbone,pbtwo,ρ,ch)
    A = pa*pbone
    return tr(pbtwo*ch*A*ρ*A*adjoint(ch))
end

function qtnumlog(pa,pbone,ρ, prob)
    return prob/(tr(pbone*pa*ρ))
end

function qtdenlog(mone,pbone,pbtwo,ρ, ch)
    A = pbtwo*ch*pbone
    B = pbone*adjoint(ch)
    a = tr(A*ρ*B)/2
    b = tr(A*mone*ρ*mone*B)/2
    c = tr(pbone*ρ)
    return (a+b)/c
end

function qtanTEsum(ma, pa,pbone,pbtwo,ρ,ch)
    prob = qtfullprob(pa, pbone, pbtwo, ρ, ch)
    if prob == 0
        return 0
    end
    return prob*log2(qtnumlog(pa, pbone, ρ, prob)/qtdenlog(ma, pbone, pbtwo, ρ, ch))
end

function qtTEan(α,βone,βtwo,ρ,ch)
    Pa = [projector(α..., i, 1) for i in [-1,1]]
    Pbone = [projector(βone..., i, 2) for i in [-1,1]]
    Pbtwo = [projector(βtwo..., i, 2) for i in [-1,1]]
    Ma = kron(spherical_point(α...),Matrix(I,2,2))
    sum = 0
    for i in Pa
        for j in Pbone
            for k in Pbtwo
                sum += real(qtanTEsum(Ma, i, j, k, ρ, ch))
            end
        end
    end
    return round(sum, digits = 4)
end

function Kraussdec(λone,λtwo,λthree)
    A = λone*kron(σx.mat,σx.mat)+λtwo*kron(σy.mat,σy.mat)+λthree*kron(σz.mat,σz.mat)
    return exp(-im*A)
end

function dmatrix(ϕ,base)
    vec = cos(ϕ)*base[1] + sin(ϕ)*base[2]
    return vec*vec'
end

function generaldm_qtTEan(α,βone,βtwo,ϕ,ch,base)
    ρ = dmatrix(ϕ,base)
    return qtTEan(α,βone,βtwo,ρ,ch)
end

function general_qtTEan(α,βone,βtwo,ϕ,λ,base)
    ch = Kraussdec(λ...)
    ρ = dmatrix(ϕ,base)
    return qtTEan(α,βone,βtwo,ρ,ch)
end

function Haar_rand_channel(n)
    A = randn(ComplexF64, n, n)
    Q,R = qr(A)
    return Q*Diagonal(R./abs.(R))
end

σx = pauli_operator([0 1; 1 0], "Z")
σy = pauli_operator([0 -im; im 0], "Z")
σz = pauli_operator([1 0; 0 -1], "Z")
Id = pauli_operator(Matrix(I, 2, 2), "Z")

qubit = [1, 2]
nqubit = length(qubit)

Mx = [measure(σx, i, nqubit) for i in qubit]
My = [measure(σy, i, nqubit) for i in qubit]
Mz = [measure(σz, i, nqubit) for i in qubit]
MId = [measure(Id, i, nqubit) for i in qubit];


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

function Dict_rot(nmeas, qubit)
     Dict("from_Y_to_X" => [measure(operator(rot(σy,σz, -(i/nmeas)*π/2), "Z"), qubit, 2) for i in 0:nmeas],
    "from_X_to_Y" => [measure(operator(rot(σx,σz, (i/nmeas)*π/2), "Z"), qubit, 2) for i in 0:nmeas],
    "from_Z_to_X" => [measure(operator(rot(σz,σy, (i/nmeas)*π/2), "Z"), qubit, 2) for i in 0:nmeas],
    "from_Z_to_Y" => [measure(operator(rot(σz,σx, -(i/nmeas)*π/2), "Z"), qubit, 2) for i in 0:nmeas],
    "from_X_to_Z" => [measure(operator(rot(σx,σy, -(i/nmeas)*π/2), "Z"), qubit, 2) for i in 0:nmeas],
    "from_Y_to_Z" => [measure(operator(rot(σy,σx, (i/nmeas)*π/2), "Z"), qubit, 2) for i in 0:nmeas]);
end