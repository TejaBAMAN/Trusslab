# julia_backend/results.jl
# Analysis result structures, summaries, and post-processing

module Results

using Printf

export NodeDisplacement, SupportReaction, ElementResult, AnalysisSummary, FullAnalysisResult,
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
end

"""
    compute_deformed_coordinates(coords::Dict{Int, Vector{Float64}}, displacements::Vector{NodeDisplacement}, scale::Float64)
"""
function compute_deformed_coordinates(coords::Dict{Int, Vector{Float64}}, displacements::Vector{NodeDisplacement}, scale::Float64)
    def_coords = Dict{Int, Vector{Float64}}()
    disp_map = Dict(d.node => (d.ux, d.uy) for d in displacements)
    
    for (node_id, pt) in coords
        ux, uy = get(disp_map, node_id, (0.0, 0.0))
        def_coords[node_id] = [pt[1] + scale * ux, pt[2] + scale * uy]
    end
    
    return def_coords
end

"""
    generate_summary(displacements::Vector{NodeDisplacement}, element_results::Vector{ElementResult}, reactions::Vector{SupportReaction}, total_fx::Float64, total_fy::Float64)
"""
function generate_summary(displacements::Vector{NodeDisplacement}, element_results::Vector{ElementResult}, reactions::Vector{SupportReaction}, total_fx::Float64, total_fy::Float64)
    # Maximum displacement
    max_disp = 0.0
    max_disp_node = 1
    for d in displacements
        if d.magnitude >= max_disp
            max_disp = d.magnitude
            max_disp_node = d.node
        end
    end
    
    # Tension / Compression / Stress
    max_tens = 0.0
    max_tens_elem = 0
    max_comp = 0.0
    max_comp_elem = 0
    max_str = 0.0
    max_str_elem = 0
    
    for el in element_results
        if el.axial_force > max_tens
            max_tens = el.axial_force
            max_tens_elem = el.id
        end
        if el.axial_force < -max_comp
            max_comp = abs(el.axial_force)
            max_comp_elem = el.id
        end
        if abs(el.stress) >= max_str
            max_str = abs(el.stress)
            max_str_elem = el.id
        end
    end
    
    # Equilibrium check
    sum_rx = sum(r.rx for r in reactions; init=0.0)
    sum_ry = sum(r.ry for r in reactions; init=0.0)
    eq_err_x = abs(total_fx + sum_rx)
    eq_err_y = abs(total_fy + sum_ry)
    eq_sat = (eq_err_x < 1e-4) && (eq_err_y < 1e-4)
    
    status = eq_sat ? "Analysis Successful: Equilibrium Satisfied" : "Analysis Finished with Minor Equilibrium Residual"
    
    return AnalysisSummary(
        max_disp, max_disp_node,
        max_tens, max_tens_elem,
        max_comp, max_comp_elem,
        max_str, max_str_elem,
        eq_err_x, eq_err_y,
        eq_sat, status
    )
end

end # module Results
