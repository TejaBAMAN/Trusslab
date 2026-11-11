# julia_backend/stability.jl
# Structural stability, determinacy, and singularity checks

module Stability

using LinearAlgebra

export StabilityReport, check_structural_stability

struct StabilityReport
    is_stable::Bool
    determinacy_status::String  # "Statically Determinate", "Statically Indeterminate", "Unstable Mechanism"
    dsi::Int                    # Degree of static indeterminacy: (b + r) - 2j
    condition_number::Float64
    num_nodes::Int
    num_elements::Int
    num_restraints::Int
    messages::Vector{String}
end

"""
    check_structural_stability(num_nodes::Int, num_elements::Int, fixed_dofs::Vector{Int}, free_dofs::Vector{Int}, K::Matrix{Float64}) -> StabilityReport

Performs static determinacy, support sufficiency, and kinematic stability checks.
"""
function check_structural_stability(num_nodes::Int, num_elements::Int, fixed_dofs::Vector{Int}, free_dofs::Vector{Int}, K::Matrix{Float64})
    j = num_nodes
    b = num_elements
    r = length(fixed_dofs)
    total_dofs = 2 * j
    dsi = (b + r) - total_dofs
    
    messages = String[]
    is_stable = true
    
    # Check 1: Minimum restraints to prevent rigid body motion in 2D plane
    if r < 3
        is_stable = false
        push!(messages, "Insufficient supports: Truss has only $r restraint(s). At least 3 non-parallel restraints are required in 2D to prevent rigid-body translations and rotation.")
    end
    
    # Check 2: Check if all restraints are in one direction (e.g. only Y or only X)
    has_x_restraint = any(d -> isodd(d), fixed_dofs)
    has_y_restraint = any(d -> iseven(d), fixed_dofs)
    if !has_x_restraint
        is_stable = false
        push!(messages, "Rigid body motion: No horizontal (X) restraint found. Truss is free to slide horizontally.")
    end
    if !has_y_restraint
        is_stable = false
        push!(messages, "Rigid body motion: No vertical (Y) restraint found. Truss is free to move vertically.")
    end
    
    # Determinacy classification
    if dsi == 0
        det_status = "Statically Determinate"
        push!(messages, "Truss is statically determinate: b + r = 2j ($b + $r = $(2*j)).")
    elseif dsi > 0
        det_status = "Statically Indeterminate"
        push!(messages, "Truss is statically indeterminate to degree $dsi: b + r > 2j ($b + $r > $(2*j)). Direct Stiffness FEM accommodates this naturally.")
    else
        det_status = "Unstable Mechanism"
        is_stable = false
        push!(messages, "Truss is statically deficient / under-constrained: b + r < 2j ($b + $r < $(2*j)). Deficit of $(abs(dsi)) member(s) or constraint(s).")
    end
    
    # Check 3: Kinematic stability via stiffness matrix condition number
    cond_num = -1.0
    if isempty(free_dofs)
        cond_num = 1.0
    elseif length(free_dofs) > 0 && is_stable
        K_ff = K[free_dofs, free_dofs]
        try
            # Compute singular values
            s_vals = svdvals(K_ff)
            min_s = s_vals[end]
            max_s = s_vals[1]
            cond_num = (min_s > 0) ? (max_s / min_s) : -1.0
            
            if min_s < 1e-8 || cond_num > 1e13
                is_stable = false
                push!(messages, "Singular stiffness matrix: Truss contains an internal mechanism or unconstrained degree of freedom (min singular value = $(round(min_s, sigdigits=4)), condition number = $(round(cond_num, sigdigits=4))).")
            else
                push!(messages, "Stiffness matrix is well-conditioned (condition number = $(round(cond_num, sigdigits=4))).")
            end
        catch err
            is_stable = false
            push!(messages, "Error factoring stiffness matrix: $(err). Structure is kinematically unstable.")
        end
    end
    
    return StabilityReport(is_stable, det_status, dsi, cond_num, j, b, r, messages)
end

end # module Stability
