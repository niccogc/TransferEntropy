#we try to generalize the tree structure to be able to use it for multiple node types

function create_tree(Circuit, density_matrix, Treetype)
    root = initializer(Treetype, Circuit, density_matrix)
    branches = length(root.children)
    for i in 1:branches
        root.children[i] = attach_node_children(root, Circuit, i)
    end
    return root
end

function attach_node_children(nod, Circuit, branch)
    if conditions(nod, Circuit) == true
        return nothing
    end
    node = operation(nod, Circuit, branch)
    branches = length(node.children)
    for i in 1:branches
        node.children[i] = attach_node_children(node, Circuit[2:end], i)
    end
    return node
end

function initializer(node::Binary_tree, Circuit, density_mat; start = 1, parent = nothing, prob = 1, res = 0)
    return Binary_tree(start, parent, density_prob(prob, density_mat.mat, res), [nothing; nothing])
end

function conditions(node::Binary_tree, Circuit)
    return length(Circuit) == 0
end

function operation(nod::Binary_tree, Circuit, branch)
    if typeof(Circuit[1]) == measure
        A = Circuit[1].proj[branch]*nod.value.mat*Circuit[1].proj[branch]
        pA = real(tr(A))
        if pA == 0
            node = Binary_tree(nod.idx*2+(branch-1), nod, density_prob(pA, zero(A), 0), [nothing; nothing])
        else
            node = Binary_tree(nod.idx*2+(branch-1), nod, density_prob(pA, A./pA, Circuit[1].value[branch]), [nothing; nothing])
        end
        return node
    elseif typeof(Circuit[1]) == channel
        return Binary_tree(nod.idx*2 +(branch-1), nod, density_prob(2-branch, Circuit[1].mat*nod.value.mat*adjoint(Circuit[1].mat), 0), [nothing; nothing])
    end
end

function trimming!(root::Binary_tree)
    if root === nothing
        return
    end
    if root.children[1] === nothing && root.children[2] === nothing
        return
    end
    if root.children[1] === nothing && root.children[2] !== nothing
        trimming!(root.children[2])
        return
    end
    if root.children[2] === nothing && root.children[1] !== nothing
        trimming!(root.children[1])
        return
    end
    if root.children[1].value.prob == 0
        root.children[1] = nothing
        trimming!(root.children[2])
        return
    end
    if root.children[2].value.prob == 0
        root.children[2] = nothing
        trimming!(root.children[1])
        return
    end
    trimming!(root.children[1])
    trimming!(root.children[2])
end

function find_node(node, target_value)
    if node === nothing
        return 
    end
    if node.idx == target_value
        return node 
    end
    A = find_node(node.children[1], target_value)
    A !== nothing ? (return A) : (return find_node(node.children[2], target_value))
end


function find_leaves(node)
    if node === nothing
        return []
    end
    if node.children[1] === nothing && node.children[2] === nothing
        return [node]
    end
    return [find_leaves(node.children[1]); find_leaves(node.children[2])]
end

function make_tree(root)
    if root === nothing
        return []
    end
    return [root; make_tree(root.children[1]); make_tree(root.children[2])]
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
        root.children[1] = node.children[1]
        root.children[2] = node.children[2]
        return
    end
    attach_node!(root.children[1], node)
    attach_node!(root.children[2], node)
end
