# julia_backend/assembly.jl
# Assembly of global stiffness matrix [K] and global load vector {F}

module Assembly

using LinearAlgebra
using ..Elements
using ..Stiffness
using ..BoundaryConditions

export LoadDef, assemble_global_stiffness, assemble_global_load_vector

struct LoadDef
    node::Int
    fx::Float64
    fy::Float64
end

"""
    assemble_global_stiffness(num_nodes::Int, elements::Vector{Element}, E::Float64, A::Float64) -> Matrix{Float64}

Assembles the (2*num_nodes) x (2*num_nodes) global stiffness matrix [K].
"""
function assemble_global_stiffness(num_nodes::Int, elements::Vector{Element}, E::Float64, A::Float64)
    total_dofs = 2 * num_nodes
    K = zeros(Float64, total_dofs, total_dofs)
    
    for el in elements
        dof_ix, dof_iy = get_node_dofs(el.node_i)
        dof_jx, dof_jy = get_node_dofs(el.node_j)
        el_dofs = [dof_ix, dof_iy, dof_jx, dof_jy]
        
        ke = calculate_element_stiffness(E, A, el.length, el.c, el.s)
        
        for a in 1:4
            for b in 1:4
                K[el_dofs[a], el_dofs[b]] += ke[a, b]
            end
        end
    end
    
    return K
end

"""
    assemble_global_load_vector(num_nodes::Int, loads::Vector{LoadDef}) -> Vector{Float64}

Assembles the (2*num_nodes) global force vector {F}.
"""
function assemble_global_load_vector(num_nodes::Int, loads::Vector{LoadDef})
    total_dofs = 2 * num_nodes
    F = zeros(Float64, total_dofs)
    
    for ld in loads
        if ld.node >= 1 && ld.node <= num_nodes
            dof_x, dof_y = get_node_dofs(ld.node)
            F[dof_x] += ld.fx
            F[dof_y] += ld.fy
        end
    end
    
    return F
end

end # module Assembly
