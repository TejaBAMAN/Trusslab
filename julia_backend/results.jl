# julia_backend/results.jl
# Analysis result structures, summaries, and post-processing
# Supports dual-axis Euler buckling evaluation per member.

module Results

using Printf

export NodeDisplacement, SupportReaction, ElementResult, ElementBucklingResult,
       AnalysisSummary, FullAnalysisResult,
       compute_deformed_coordinates, generate_summary

struct NodeDisplacement
    node::Int
    ux::Float64
    uy::Float64
    magnitude::Float64
end

struct SupportReaction
    node::Int
    rx::Float64
    ry::Float64
    magnitude::Float64
end

struct ElementResult
    id::Int
    node_i::Int
    node_j::Int
    length::Float64
    c::Float64
    s::Float64
    angle_deg::Float64
    elongation::Float64
    strain::Float64
    stress::Float64
    axial_force::Float64
    state::String # "Tension", "Compression", "Zero Force"
end

struct ElementBucklingResult
    id::Int
    section_name::String
    section_shape::String
    A::Float64            # m^2
    Iu::Float64           # m^4
    Iv::Float64           # m^4
    ku::Float64
    kv::Float64
    Lu::Float64           # m
    Lv::Float64           # m
    p_cr_u::Float64       # N
    p_cr_v::Float64       # N
    p_cr_gov::Float64     # N
    governing_axis::String # "u", "v", or "N/A"
    axial_force::Float64  # N
    utilisation::Float64  # |N| / Pcr_gov
    status::String        # "Safe", "Buckled", "Tension (N/A)", "Zero Force (N/A)", "Not Configured"
    configured::Bool
end

struct AnalysisSummary
    max_displacement::Float64
    max_displacement_node::Int
    max_tension::Float64
    max_tension_element::Int
    max_compression::Float64
    max_compression_element::Int
    max_stress::Float64
    max_stress_element::Int
    equilibrium_error_x::Float64
    equilibrium_error_y::Float64
    equilibrium_satisfied::Bool
    status::String
end

struct FullAnalysisResult
    success::Bool
    error_message::String
    stability_report::Any
    displacements::Vector{NodeDisplacement}
    reactions::Vector{SupportReaction}
    elements::Vector{ElementResult}
    summary::AnalysisSummary
    deformed_coords::Dict{Int, Vector{Float64}}
    buckling_results::Vector{ElementBucklingResult}
end

# Backward compatibility constructor for 8-argument callers
function FullAnalysisResult(
    success::Bool,
    error_message::String,
    stability_report::Any,
    displacements::Vector{NodeDisplacement},
    reactions::Vector{SupportReaction},
    elements::Vector{ElementResult},
    summary::AnalysisSummary,
    deformed_coords::Dict{Int, Vector{Float64}}
)
    return FullAnalysisResult(
        success, error_message, stability_report,
        displacements, reactions, elements, summary, deformed_coords,
        ElementBucklingResult[]
    )
end

"""
    compute_deformed_coordinates(coords::Dict{Int, Vector{Float64}}, displacements::Vector{NodeDisplacement}, scale::Float64)
"""
function compute_deformed_coordinates(coords::Dict{Int, Vector{Float64}}, displacements::Vector{NodeDisplacement}, scale::Float64)
    def_coords = Dict{Int, Vector{Float64}}()
    disp_map = Dict(d.node => [d.ux, d.uy] for d in displacements)
    
    for (node, pt) in coords
        d = get(disp_map, node, [0.0, 0.0])
        def_coords[node] = [pt[1] + scale * d[1], pt[2] + scale * d[2]]
    end
    
    return def_coords
end

"""
    generate_summary(displacements, elements, reactions, total_fx, total_fy; tol=1e-4) -> AnalysisSummary
"""
function generate_summary(
    displacements::Vector{NodeDisplacement},
    elements::Vector{ElementResult},
    reactions::Vector{SupportReaction},
    total_fx::Float64,
    total_fy::Float64;
    tol::Float64=1e-4
)
    max_disp = 0.0
    max_disp_node = 0
    for d in displacements
        if d.magnitude > max_disp
            max_disp = d.magnitude
            max_disp_node = d.node
        end
    end
    
    max_tens = 0.0
    max_tens_elem = 0
    max_comp = 0.0
    max_comp_elem = 0
    max_stress = 0.0
    max_stress_elem = 0
    
    for el in elements
        if el.axial_force > max_tens
            max_tens = el.axial_force
            max_tens_elem = el.id
        end
        if el.axial_force < -max_comp
            max_comp = abs(el.axial_force)
            max_comp_elem = el.id
        end
        abs_s = abs(el.stress)
        if abs_s > max_stress
            max_stress = abs_s
            max_stress_elem = el.id
        end
    end
    
    sum_rx = sum(r.rx for r in reactions; init=0.0)
    sum_ry = sum(r.ry for r in reactions; init=0.0)
    
    err_x = abs(total_fx + sum_rx)
    err_y = abs(total_fy + sum_ry)
    eq_satisfied = (err_x <= tol) && (err_y <= tol)
    
    status = eq_satisfied ? "Converged & In Equilibrium" : "Equilibrium Not Satisfied"
    
    return AnalysisSummary(
        max_disp, max_disp_node,
        max_tens, max_tens_elem,
        max_comp, max_comp_elem,
        max_stress, max_stress_elem,
        err_x, err_y, eq_satisfied, status
    )
end

end # module Results
