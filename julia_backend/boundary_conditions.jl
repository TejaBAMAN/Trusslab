# julia_backend/boundary_conditions.jl
# Support conditions and DOF partitioning

module BoundaryConditions

export SupportDef, identify_dofs, get_node_dofs

struct SupportDef
    node::Int
    type::String  # "pinned", "roller_x" (UY fixed), "roller_y" (UX fixed), "fixed", "custom"
    fix_x::Bool
    fix_y::Bool
end

"""
    get_node_dofs(node::Int) -> (Int, Int)
Returns global DOF indices (dof_x, dof_y) for a node.
"""
function get_node_dofs(node::Int)
    return (2 * node - 1, 2 * node)
end

"""
    identify_dofs(num_nodes::Int, supports::Vector{SupportDef}) -> (fixed_dofs::Vector{Int}, free_dofs::Vector{Int})
"""
function identify_dofs(num_nodes::Int, supports::Vector{SupportDef})
    fixed_set = Set{Int}()
    
    for s in supports
        dof_x, dof_y = get_node_dofs(s.node)
        if s.fix_x
            push!(fixed_set, dof_x)
        end
        if s.fix_y
            push!(fixed_set, dof_y)
        end
    end
    
    total_dofs = 2 * num_nodes
    fixed_dofs = sort(collect(fixed_set))
    free_dofs = [d for d in 1:total_dofs if !(d in fixed_set)]
    
    return (fixed_dofs, free_dofs)
end

end # module BoundaryConditions
