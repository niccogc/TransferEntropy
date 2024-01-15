
function children(node::Binary_tree)
    if node === nothing
        return []
    end
    return [node.left; node.right], [1,2]
end

function attach_node_left(meas::Vector{operators}, idx, mat, nod)
    length(meas) == 0 && return nothing
    #println(typeof(meas[1].proj[1]))
    #println(typeof(mat))
    node = node_operator(meas[1], idx, mat, nod, 1)
    node.left = attach_node_left(meas[2:end], node.idx, node.value.mat, node)
    node.right = attach_node_right(meas[2:end], node.idx, node.value.mat, node)
    return node
end

function attach_node_right(meas::Vector{operators}, idx, mat, nod)
    length(meas) == 0 && return nothing
    node = node_operator(meas[1], idx, mat, nod, 2)
    node.left = attach_node_left(meas[2:end], node.idx, node.value.mat, node)
    node.right = attach_node_right(meas[2:end], node.idx, node.value.mat, node)
    return node
end

function node_operator(measure::measure, idx, mat, nod, branch)
    A = measure.proj[branch]*mat*measure.proj[branch]
    #println(A)
    pA = tr(A)
    if pA == 0
        node = Binary_tree(idx*2+(branch-1), nod, density_prob(pA, zero(A), 0), nothing, nothing)
    else
        node = Binary_tree(idx*2+(branch-1), nod, density_prob(pA, A./pA, measure.value[branch]), nothing, nothing)
    end
    return node
end

function node_operator(chan::channel, idx, mat, nod, branch)
    return Binary_tree(idx*2 +(branch-1), nod, density_prob(branch-1, chan.mat*mat*adjoint(chan.mat), 0), nothing, nothing)
end

function create_tree(meas::Vector{operators}, density_mat::density_matrix; start = 1, parent = nothing, prob = 1, res = 0)
    root = Binary_tree(start, parent, density_prob(prob, density_mat.mat, res), nothing, nothing)
    root.left = attach_node_left(meas, root.idx, root.value.mat, root)
    root.right = attach_node_right(meas, root.idx, root.value.mat, root)
    return root
end


function trimming!(root)
    if root === nothing
        return
    end
    #=
    if root.value.prob == 0
        root.value = density_prob(0, zero(root.value.mat), 0)
        root.left = nothing
        root.right = nothing
        return
    end
    =#
    if root.left === nothing && root.right === nothing
        return
    end
    if root.left === nothing && root.right !== nothing
        trimming!(root.right)
        return
    end
    if root.right === nothing && root.left !== nothing
        trimming!(root.left)
        return
    end
    if root.left.value.prob == 0
        root.left = nothing
        trimming!(root.right)
        return
    end
    if root.right.value.prob == 0
        root.right = nothing
        trimming!(root.left)
        return
    end
    trimming!(root.left)
    trimming!(root.right)
end

function find_node(node, target_value)
    if node === nothing
        return 
    end
    if node.idx == target_value
        return node 
    end
    A = find_node(node.left, target_value)
    A !== nothing ? (return A) : (return find_node(node.right, target_value))
end


function find_leaves(node)
    if node === nothing
        return []
    end
    if node.left === nothing && node.right === nothing
        return [node]
    end
    return [find_leaves(node.left); find_leaves(node.right)]
end

#=
function make_tree(root::TreeNode)
    A = []
    i=0
    while true
        i += 1
        J = find_node(root, i)
        if J === nothing
            push!(A, J)
        else
            push!(A,J)
            if J.left === nothing && J.right === nothing
                break
            end
        end
    end
    return Tree(A)
end
=#
function make_tree(root)
    if root === nothing
        return []
    end
    return [root; make_tree(root.left); make_tree(root.right)]
end

function find_root_path(node::TreeNode)
    if node === nothing
        return []
    end
    if node.parent === nothing
        return [node]
    end
    return [node; find_root_path(node.parent)]
end

function attach_node!(root, node)
    if root === nothing
        return
    end
    if root.idx == node.idx
        root.idx = node.idx
        root.parent = node.parent
        root.value = node.value
        root.left = node.left
        root.right = node.right
        return
    end
    attach_node!(root.left, node)
    attach_node!(root.right, node)
end