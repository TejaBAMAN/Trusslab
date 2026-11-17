# julia_backend/solver.jl
# Finite Element Direct Stiffness solver for 2D trusses
# Supports heterogeneous per-member section properties and dual-axis Euler buckling verification.

module Solver

using LinearAlgebra
using Printf
using ..Elements
using ..Stiffness
using ..BoundaryConditions
using ..Assembly
using ..Stability
using ..Results

export solve_truss_fem

"""
    solve_truss_fem(num_nodes::Int, elements::Vector{Element}, coords::Dict{Int, Vector{Float64}}, E::Float64, A::Float64, supports::Vector{SupportDef}, loads::Vector{LoadDef}; ...) -> FullAnalysisResult

Solves the 2D truss using the Finite Element Direct Stiffness Method.
Supports:
- Formulations: "dsm", "cholesky", "lu", "qr"
- Analysis Types: "linear_elastic", "second_order", "buckling", "gravity"
- Heterogeneous per-member properties via `member_configs` parameter
- Dual-axis Euler elastic buckling verification about principal u and v axes
"""
function solve_truss_fem(
    num_nodes::Int,
    elements::Vector{Element},
    coords::Dict{Int, Vector{Float64}},
    E::Float64,
    A::Float64,
    supports::Vector{SupportDef},
    loads::Vector{LoadDef};
    member_configs::Vector{MemberConfig}=MemberConfig[],
    def_scale::Float64=100.0,
    analysis_type::String="linear_elastic",
    formulation::String="dsm",
    density::Float64=7850.0,
    include_self_weight::Bool=false,
    check_buckling::Bool=false
)
    # Determine flags with backward compatibility
    use_gravity = include_self_weight || (analysis_type == "gravity")
    eval_buckling = check_buckling || (analysis_type == "buckling")
    base_analysis = (analysis_type == "second_order") ? "second_order" : "linear_elastic"

    num_elems = length(elements)
    
    # Prepare per-member properties
    cfgs = MemberConfig[]
    if length(member_configs) == num_elems
        cfgs = member_configs
    else
        # Fall back to uniform scalar E, A, density with default circular section and pinned-pinned buckling
        default_sec = compute_solid_circle(sqrt(4.0 * A / pi); name="Default Circular Bar")
        for el in elements
            push!(cfgs, MemberConfig(el.id, default_sec, BucklingRestraint("pinned_pinned", "pinned_pinned", 1.0, 1.0, el.length, el.length, eval_buckling), E, density))
        end
    end

    Es = [c.E for c in cfgs]
    As = [c.section.A for c in cfgs]
    densities = [c.density for c in cfgs]

    # Step 1: Assemble Global Stiffness Matrix and Force Vector
    K = assemble_global_stiffness(num_nodes, elements, Es, As)
    F = assemble_global_load_vector(num_nodes, loads)
    
    # Gravitational Dead Load inclusion if requested (using individual member areas & densities)
    total_dead_load = 0.0
    if use_gravity
        g = 9.80665 # m/s^2
        for (idx, el) in enumerate(elements)
            if densities[idx] > 0.0 && As[idx] > 0.0
                el_mass = densities[idx] * As[idx] * el.length
                el_weight = el_mass * g
                total_dead_load += el_weight
                dof_iy = 2 * el.node_i
                dof_jy = 2 * el.node_j
                F[dof_iy] -= 0.5 * el_weight
                F[dof_jy] -= 0.5 * el_weight
            end
        end
    end
    
    # Step 2: Partition Boundary Conditions
    fixed_dofs, free_dofs = identify_dofs(num_nodes, supports)
    
    # Step 3: Stability & Singularity Check
    stab_report = check_structural_stability(num_nodes, length(elements), fixed_dofs, free_dofs, K)
    if !stab_report.is_stable
        return FullAnalysisResult(
            false,
            join(stab_report.messages, " \n"),
            stab_report,
            NodeDisplacement[],
            SupportReaction[],
            ElementResult[],
            AnalysisSummary(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, false, "Failed Stability Check"),
            coords,
            ElementBucklingResult[]
        )
    end
    
    # Step 4: Solve for Free Displacements: [K_ff] {U_f} = {F_f}
    total_dofs = 2 * num_nodes
    U = zeros(Float64, total_dofs)
    
    if !isempty(free_dofs)
        K_ff = K[free_dofs, free_dofs]
        F_f = F[free_dofs]
        try
            norm_form = lowercase(strip(formulation))
            if norm_form == "cholesky"
                # Cholesky LL^T decomposition for Symmetric Positive Definite stiffness
                K_sym = Symmetric(K_ff)
                U_f = cholesky(K_sym) \ F_f
            elseif norm_form == "lu"
                # LU factorization with partial pivoting
                U_f = lu(K_ff) \ F_f
            elseif norm_form == "qr"
                # QR orthogonal-triangular factorization
                U_f = qr(K_ff) \ F_f
            else
                # Default Direct Stiffness (DSM) linear solve
                U_f = K_ff \ F_f
            end

            if any(isnan.(U_f)) || any(isinf.(U_f))
                error("Displacements are infinite or NaN (Singular/Ill-conditioned stiffness matrix).")
            end
            U[free_dofs] = U_f
        catch err
            return FullAnalysisResult(
                false,
                "Numerical error solving [K]{U} = {F} via $(uppercase(formulation)): $(err)",
                stab_report,
                NodeDisplacement[],
                SupportReaction[],
                ElementResult[],
                AnalysisSummary(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, false, "Solver Factorization Error"),
                coords,
                ElementBucklingResult[]
            )
        end
    end
    
    # Step 4b: Geometric Non-linear (P-Delta) Iteration if requested
    if base_analysis == "second_order" && !isempty(free_dofs)
        try
            coords_def = Dict{Int, Vector{Float64}}()
            for n in 1:num_nodes
                coords_def[n] = [coords[n][1] + U[2*n - 1], coords[n][2] + U[2*n]]
            end
            elems_def = Element[]
            for el in elements
                p1 = coords_def[el.node_i]
                p2 = coords_def[el.node_j]
                dx = p2[1] - p1[1]
                dy = p2[2] - p1[2]
                L = hypot(dx, dy)
                push!(elems_def, Element(el.id, el.node_i, el.node_j, L, dx/L, dy/L, atan(dy, dx), rad2deg(atan(dy, dx))))
            end
            K2 = assemble_global_stiffness(num_nodes, elems_def, Es, As)
            K2_ff = K2[free_dofs, free_dofs]
            U2_f = K2_ff \ F[free_dofs]
            if !any(isnan.(U2_f)) && !any(isinf.(U2_f))
                U[free_dofs] = U2_f
                K = K2
            end
        catch
            # Retain first-order displacements if second-order update fails
        end
    end

    # Step 5: Calculate Support Reactions: {R} = [K]{U} - {F}
    R_total = K * U - F
    
    # Extract Nodal Displacements
    displacements = NodeDisplacement[]
    for n in 1:num_nodes
        dof_x, dof_y = get_node_dofs(n)
        ux = U[dof_x]
        uy = U[dof_y]
        push!(displacements, NodeDisplacement(n, ux, uy, hypot(ux, uy)))
    end
    
    # Extract Support Reactions (at restrained nodes)
    reactions = SupportReaction[]
    for s in supports
        dof_x, dof_y = get_node_dofs(s.node)
        rx = s.fix_x ? R_total[dof_x] : 0.0
        ry = s.fix_y ? R_total[dof_y] : 0.0
        push!(reactions, SupportReaction(s.node, rx, ry, hypot(rx, ry)))
    end
    
    # Step 6: Calculate Element Elongation, Strain, Stress, Axial Force and Dual-Axis Buckling
    element_results = ElementResult[]
    buckling_results = ElementBucklingResult[]
    
    for (idx, el) in enumerate(elements)
        ni = el.node_i
        nj = el.node_j
        
        dof_ix, dof_iy = get_node_dofs(ni)
        dof_jx, dof_jy = get_node_dofs(nj)
        
        uix = U[dof_ix]
        uiy = U[dof_iy]
        ujx = U[dof_jx]
        ujy = U[dof_jy]
        
        c = el.c
        s = el.s
        L = el.length
        
        # Elongation: Delta L = c * (u_jx - u_ix) + s * (u_jy - u_iy)
        delta_L = c * (ujx - uix) + s * (ujy - uiy)
        strain = delta_L / L
        axial_force = delta_L * (Es[idx] * As[idx] / L)
        stress = axial_force / As[idx]
        
        # Dual-Axis Euler Buckling Verification
        cfg = cfgs[idx]
        sec = cfg.section
        bck = cfg.buckling
        
        Iu = sec.Iu
        Iv = sec.Iv
        ku = bck.ku > 0.0 ? bck.ku : 1.0
        kv = bck.kv > 0.0 ? bck.kv : 1.0
        Lu = bck.Lu > 0.0 ? bck.Lu : L
        Lv = bck.Lv > 0.0 ? bck.Lv : L
        
        p_cr_u = (pi^2 * Es[idx] * Iu) / ((ku * Lu)^2)
        p_cr_v = (pi^2 * Es[idx] * Iv) / ((kv * Lv)^2)
        p_cr_gov = min(p_cr_u, p_cr_v)
        gov_axis = (p_cr_u <= p_cr_v) ? "u" : "v"
        
        util = 0.0
        status = "Not Configured"
        state_str = "Zero Force"
        
        if axial_force > 1e-6
            state_str = "Tension"
            status = "Tension (N/A)"
            util = 0.0
        elseif axial_force < -1e-6
            if eval_buckling && bck.configured
                util = abs(axial_force) / p_cr_gov
                if util >= 1.0
                    status = "Buckled"
                    state_str = @sprintf("Buckled (%.0f%% Pcr - %s)", util * 100, gov_axis)
                else
                    status = "Safe"
                    state_str = @sprintf("Compression (%.0f%% Pcr - %s)", util * 100, gov_axis)
                end
            else
                status = bck.configured ? "Safe" : "Not Configured"
                state_str = "Compression"
            end
        else
            state_str = "Zero Force"
            status = "Zero Force (N/A)"
            util = 0.0
        end
        
        push!(element_results, ElementResult(
            el.id, ni, nj, L, c, s, el.angle_deg, delta_L, strain, stress, axial_force, state_str
        ))
        
        push!(buckling_results, ElementBucklingResult(
            el.id, sec.name, sec.shape, As[idx], Iu, Iv, ku, kv, Lu, Lv,
            p_cr_u, p_cr_v, p_cr_gov, gov_axis, axial_force, util, status, bck.configured
        ))
    end
    
    # Step 7: Global Equilibrium & Summary
    total_fx = sum(ld.fx for ld in loads; init=0.0)
    total_fy = sum(ld.fy for ld in loads; init=0.0) - total_dead_load
    summary = generate_summary(displacements, element_results, reactions, total_fx, total_fy)
    
    # Status augmentation
    status_parts = String[]
    if base_analysis == "second_order"
        push!(status_parts, "Geometric Nonlinear (P-Delta)")
    else
        push!(status_parts, "Linear Elastic (DSM)")
    end
    if use_gravity
        push!(status_parts, @sprintf("+ Dead Load (%.2f kN)", total_dead_load / 1000.0))
    end
    if eval_buckling
        configured_comp = filter(b -> b.axial_force < -1e-6 && b.configured, buckling_results)
        if !isempty(configured_comp)
            max_util = maximum([b.utilisation for b in configured_comp])
            lambda_cr = max_util > 0 ? (1.0 / max_util) : Inf
            push!(status_parts, @sprintf("• Dual-Axis Buckling (lambda_cr = %.2f×)", lambda_cr))
        else
            push!(status_parts, "• Buckling Checked (No active compression members)")
        end
    end
    
    summary = AnalysisSummary(
        summary.max_displacement, summary.max_displacement_node,
        summary.max_tension, summary.max_tension_element,
        summary.max_compression, summary.max_compression_element,
        summary.max_stress, summary.max_stress_element,
        summary.equilibrium_error_x, summary.equilibrium_error_y,
        summary.equilibrium_satisfied,
        join(status_parts, " ")
    )
    
    deformed_coords = compute_deformed_coordinates(coords, displacements, def_scale)
    
    return FullAnalysisResult(
        true,
        "",
        stab_report,
        displacements,
        reactions,
        element_results,
        summary,
        deformed_coords,
        buckling_results
    )
end

end # module Solver
