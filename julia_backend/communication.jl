# julia_backend/communication.jl
# JSON input/output communication layer for GUI and CLI integration
# Supports per-member cross-section and dual-axis Euler buckling payloads.

module Communication

using JSON
using ..Validation
using ..Geometry
using ..Triangles
using ..Elements
using ..Stiffness
using ..BoundaryConditions
using ..Assembly
using ..Stability
using ..Results
using ..Solver

export handle_request_json, parse_and_process

function dict_to_triangles(raw_triangles::Vector{<:Any})
    triangles = TriangleDef[]
    for (i, t) in enumerate(raw_triangles)
        tid = get(t, "id", i)
        nodes = Int.(t["nodes"])
        lens = Float64.(t["lengths"])
        side = get(t, "side", "auto")
        push!(triangles, TriangleDef(tid, nodes, lens, side))
    end
    return triangles
end

function dict_to_supports(raw_supports::Vector{<:Any})
    supports = SupportDef[]
    for s in raw_supports
        node = Int(s["node"])
        stype = get(s, "type", "custom")
        fix_x = get(s, "fix_x", false)
        fix_y = get(s, "fix_y", false)
        if stype == "pinned" || stype == "fixed"
            fix_x = true
            fix_y = true
        elseif stype == "roller_x" # Y fixed, free in X
            fix_x = false
            fix_y = true
        elseif stype == "roller_y" # X fixed, free in Y
            fix_x = true
            fix_y = false
        end
        push!(supports, SupportDef(node, stype, fix_x, fix_y))
    end
    return supports
end

function dict_to_loads(raw_loads::Vector{<:Any})
    loads = LoadDef[]
    for ld in raw_loads
        node = Int(ld["node"])
        fx = Float64(get(ld, "fx", 0.0))
        fy = Float64(get(ld, "fy", 0.0))
        push!(loads, LoadDef(node, fx, fy))
    end
    return loads
end

"""
    handle_request_json(json_str::String) -> String
Processes JSON request string and returns JSON response string.
"""
function handle_request_json(json_str::String)
    req = JSON.parse(json_str)
    res = parse_and_process(req)
    return JSON.json(res)
end

function parse_and_process(req::AbstractDict)
    action = get(req, "action", "full_pipeline")
    
    try
        if action == "validate_triangle"
            a = Float64(req["a"])
            b = Float64(req["b"])
            c = Float64(req["c"])
            val, msg = validate_triangle_lengths(a, b, c)
            angles = val ? calculate_triangle_angles(a, b, c) : (0.0, 0.0, 0.0)
            return Dict{String, Any}(
                "success" => val,
                "message" => msg,
                "angles_deg" => [rad2deg(angles[1]), rad2deg(angles[2]), rad2deg(angles[3])]
            )
            
        elseif action == "reconstruct_geometry"
            num_nodes = Int(req["num_nodes"])
            raw_tris = req["triangles"]
            
            # Validate connectivity first
            val_conn, msg_conn = validate_connectivity(num_nodes, raw_tris)
            if !val_conn
                return Dict{String, Any}("success" => false, "error" => msg_conn)
            end
            
            triangles = dict_to_triangles(raw_tris)
            coords, history = reconstruct_geometry_from_triangles(num_nodes, triangles)
            elements = generate_unique_elements(triangles, coords)
            
            elem_list = [
                Dict{String, Any}(
                    "id" => el.id,
                    "node_i" => el.node_i,
                    "node_j" => el.node_j,
                    "length" => el.length,
                    "c" => el.c,
                    "s" => el.s,
                    "angle_rad" => el.angle_rad,
                    "angle_deg" => el.angle_deg
                ) for el in elements
            ]
            
            coords_dict = Dict{String, Vector{Float64}}(string(k) => v for (k, v) in coords)
            
            return Dict{String, Any}(
                "success" => true,
                "num_nodes" => num_nodes,
                "num_triangles" => length(triangles),
                "num_elements" => length(elements),
                "coordinates" => coords_dict,
                "elements" => elem_list,
                "triangle_history" => history
            )
            
        elseif action == "analyze" || action == "full_pipeline"
            num_nodes = Int(req["num_nodes"])
            raw_tris = req["triangles"]
            
            # 1. Geometry reconstruction
            val_conn, msg_conn = validate_connectivity(num_nodes, raw_tris)
            if !val_conn
                return Dict{String, Any}("success" => false, "error" => msg_conn)
            end
            
            triangles = dict_to_triangles(raw_tris)
            coords, history = reconstruct_geometry_from_triangles(num_nodes, triangles)
            elements = generate_unique_elements(triangles, coords)
            
            # 2. Material properties
            mat = req["material"]
            E = Float64(mat["E"])
            A = Float64(mat["A"])
            density = Float64(get(mat, "density", 7850.0))
            
            # 3. Supports & Loads
            supports = dict_to_supports(req["supports"])
            loads = dict_to_loads(req["loads"])
            def_scale = Float64(get(req, "def_scale", 100.0))
            analysis_type = String(get(req, "analysis_type", "linear_elastic"))
            formulation = String(get(req, "formulation", "dsm"))
            include_self_weight = Bool(get(req, "include_self_weight", false))
            check_buckling = Bool(get(req, "check_buckling", false))
            
            # 4. Optional Per-Member Section & Buckling Configurations
            member_configs = MemberConfig[]
            raw_members = get(req, "member_properties", get(req, "member_configs", Any[]))
            
            if !isempty(raw_members)
                for (idx, el) in enumerate(elements)
                    # Find matching configuration by ID or index
                    m_data = nothing
                    for rm in raw_members
                        if get(rm, "id", get(rm, "element_id", -1)) == el.id
                            m_data = rm
                            break
                        end
                    end
                    if m_data === nothing && idx <= length(raw_members)
                        m_data = raw_members[idx]
                    end
                    
                    if m_data !== nothing
                        sec_dict = get(m_data, "section", Dict{String, Any}())
                        shape = get(sec_dict, "shape", "solid_circle")
                        sec_name = get(sec_dict, "name", "Member $(el.id) Section")
                        
                        raw_params = get(sec_dict, "params", Dict{String, Any}())
                        params = Dict{String, Float64}()
                        for (pk, pv) in raw_params
                            params[string(pk)] = Float64(pv)
                        end
                        if !haskey(params, "A") && haskey(m_data, "A")
                            params["A"] = Float64(m_data["A"])
                        end
                        if isempty(params)
                            params["D"] = sqrt(4.0 * A / pi)
                        end
                        
                        sec = compute_section_properties(shape, params; name=sec_name)
                        
                        bck_dict = get(m_data, "buckling", Dict{String, Any}())
                        preset_u = String(get(bck_dict, "preset_u", "pinned_pinned"))
                        preset_v = String(get(bck_dict, "preset_v", "pinned_pinned"))
                        ku = Float64(get(bck_dict, "ku", 1.0))
                        kv = Float64(get(bck_dict, "kv", 1.0))
                        Lu = Float64(get(bck_dict, "Lu", el.length))
                        Lv = Float64(get(bck_dict, "Lv", el.length))
                        configured = Bool(get(bck_dict, "configured", check_buckling))
                        
                        bck = BucklingRestraint(preset_u, preset_v, ku, kv, Lu, Lv, configured)
                        m_E = Float64(get(m_data, "E", E))
                        m_density = Float64(get(m_data, "density", density))
                        
                        push!(member_configs, MemberConfig(el.id, sec, bck, m_E, m_density))
                    else
                        # Default fallback for this element
                        default_sec = compute_solid_circle(sqrt(4.0 * A / pi); name="Default Section")
                        default_bck = BucklingRestraint("pinned_pinned", "pinned_pinned", 1.0, 1.0, el.length, el.length, check_buckling)
                        push!(member_configs, MemberConfig(el.id, default_sec, default_bck, E, density))
                    end
                end
            end
            
            # 5. Solve FEM
            fem_result = solve_truss_fem(
                num_nodes, elements, coords, E, A, supports, loads;
                member_configs=member_configs,
                def_scale=def_scale,
                analysis_type=analysis_type,
                formulation=formulation,
                density=density,
                include_self_weight=include_self_weight,
                check_buckling=check_buckling
            )
            
            if !fem_result.success
                return Dict{String, Any}(
                    "success" => false,
                    "error" => fem_result.error_message,
                    "stability_report" => Dict{String, Any}(
                        "is_stable" => fem_result.stability_report.is_stable,
                        "status" => fem_result.stability_report.determinacy_status,
                        "dsi" => fem_result.stability_report.dsi,
                        "condition_number" => fem_result.stability_report.condition_number,
                        "messages" => fem_result.stability_report.messages
                    )
                )
            end
            
            # 6. Package results
            coords_dict = Dict{String, Vector{Float64}}(string(k) => v for (k, v) in coords)
            def_coords_dict = Dict{String, Vector{Float64}}(string(k) => v for (k, v) in fem_result.deformed_coords)
            
            elem_list = [
                Dict{String, Any}(
                    "id" => el.id,
                    "node_i" => el.node_i,
                    "node_j" => el.node_j,
                    "length" => el.length,
                    "c" => el.c,
                    "s" => el.s,
                    "angle_deg" => el.angle_deg,
                    "elongation" => el.elongation,
                    "strain" => el.strain,
                    "stress" => el.stress,
                    "axial_force" => el.axial_force,
                    "state" => el.state
                ) for el in fem_result.elements
            ]
            
            buckling_list = [
                Dict{String, Any}(
                    "id" => b.id,
                    "section_name" => b.section_name,
                    "section_shape" => b.section_shape,
                    "A" => b.A,
                    "Iu" => b.Iu,
                    "Iv" => b.Iv,
                    "ku" => b.ku,
                    "kv" => b.kv,
                    "Lu" => b.Lu,
                    "Lv" => b.Lv,
                    "p_cr_u" => b.p_cr_u,
                    "p_cr_v" => b.p_cr_v,
                    "p_cr_gov" => b.p_cr_gov,
                    "governing_axis" => b.governing_axis,
                    "axial_force" => b.axial_force,
                    "utilisation" => b.utilisation,
                    "status" => b.status,
                    "configured" => b.configured
                ) for b in fem_result.buckling_results
            ]
            
            disp_list = [
                Dict{String, Any}(
                    "node" => d.node,
                    "ux" => d.ux,
                    "uy" => d.uy,
                    "magnitude" => d.magnitude
                ) for d in fem_result.displacements
            ]
            
            react_list = [
                Dict{String, Any}(
                    "node" => r.node,
                    "rx" => r.rx,
                    "ry" => r.ry,
                    "magnitude" => r.magnitude
                ) for r in fem_result.reactions
            ]
            
            summary_dict = Dict{String, Any}(
                "max_displacement" => fem_result.summary.max_displacement,
                "max_displacement_node" => fem_result.summary.max_displacement_node,
                "max_tension" => fem_result.summary.max_tension,
                "max_tension_element" => fem_result.summary.max_tension_element,
                "max_compression" => fem_result.summary.max_compression,
                "max_compression_element" => fem_result.summary.max_compression_element,
                "max_stress" => fem_result.summary.max_stress,
                "max_stress_element" => fem_result.summary.max_stress_element,
                "equilibrium_error_x" => fem_result.summary.equilibrium_error_x,
                "equilibrium_error_y" => fem_result.summary.equilibrium_error_y,
                "equilibrium_satisfied" => fem_result.summary.equilibrium_satisfied,
                "status" => fem_result.summary.status
            )
            
            stab_dict = Dict{String, Any}(
                "is_stable" => fem_result.stability_report.is_stable,
                "determinacy_status" => fem_result.stability_report.determinacy_status,
                "dsi" => fem_result.stability_report.dsi,
                "condition_number" => fem_result.stability_report.condition_number,
                "messages" => fem_result.stability_report.messages
            )
            
            return Dict{String, Any}(
                "success" => true,
                "analysis_type" => analysis_type,
                "formulation" => formulation,
                "include_self_weight" => include_self_weight,
                "check_buckling" => check_buckling,
                "num_nodes" => num_nodes,
                "num_triangles" => length(triangles),
                "num_elements" => length(elements),
                "coordinates" => coords_dict,
                "deformed_coordinates" => def_coords_dict,
                "elements" => elem_list,
                "buckling_results" => buckling_list,
                "displacements" => disp_list,
                "reactions" => react_list,
                "summary" => summary_dict,
                "stability" => stab_dict
            )
        else
            return Dict{String, Any}("success" => false, "error" => "Unknown action: $action")
        end
    catch err
        return Dict{String, Any}(
            "success" => false,
            "error" => "Julia Engine Exception: $(err)",
            "stacktrace" => sprint(showerror, err, catch_backtrace())
        )
    end
end

end # module Communication
