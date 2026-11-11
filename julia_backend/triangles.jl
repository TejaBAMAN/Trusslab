# julia_backend/triangles.jl
# Progressive triangle construction, shared-edge detection, and coordinate reconstruction

module Triangles

using ..Geometry
using ..Validation

export TriangleDef, reconstruct_geometry_from_triangles, find_shared_edge_info

struct TriangleDef
    id::Int
    nodes::Vector{Int}
    lengths::Vector{Float64} # [L12, L23, L31]
    side::String             # "auto", "left", "right"
end

"""
    find_shared_edge_info(tri_nodes::Vector{Int}, tri_lengths::Vector{Float64}, coords::Dict{Int, Vector{Float64}})

Finds the shared edge between the new triangle and existing placed nodes.
Returns:
- (A::Int, B::Int): start and end nodes of directed shared edge
- C::Int: the new unplaced node
- c_len::Float64: length of shared edge AB
- b_len::Float64: length of edge AC
- a_len::Float64: length of edge BC
"""
function find_shared_edge_info(tri_nodes::Vector{Int}, tri_lengths::Vector{Float64}, coords::Dict{Int, Vector{Float64}})
    N1, N2, N3 = tri_nodes[1], tri_nodes[2], tri_nodes[3]
    L12, L23, L31 = tri_lengths[1], tri_lengths[2], tri_lengths[3]
    
    known = [n for n in tri_nodes if haskey(coords, n)]
    
    if length(known) == 2
        # Exactly one edge is shared
        if !haskey(coords, N3)
            # N1 and N2 are known, N3 is new
            A = N1
            B = N2
            C = N3
            c_len = L12
            b_len = L31 # AC
            a_len = L23 # BC
            return (A, B, C, c_len, b_len, a_len)
        elseif !haskey(coords, N1)
            # N2 and N3 are known, N1 is new
            A = N2
            B = N3
            C = N1
            c_len = L23
            b_len = L12 # AC (N2-N1)
            a_len = L31 # BC (N3-N1)
            return (A, B, C, c_len, b_len, a_len)
        else
            # N3 and N1 are known, N2 is new
            A = N3
            B = N1
            C = N2
            c_len = L31
            b_len = L23 # AC (N3-N2)
            a_len = L12 # BC (N1-N2)
            return (A, B, C, c_len, b_len, a_len)
        end
    elseif length(known) == 3
        # All three nodes already exist (closure triangle / braced panel)
        return (N1, N2, N3, L12, L31, L23)
    else
        error("Cannot construct triangle: at least 2 nodes must already have coordinates or be part of Triangle 1.")
    end
end

"""
    reconstruct_geometry_from_triangles(num_nodes::Int, triangles::Vector{TriangleDef})

Progressively reconstructs all node coordinates from triangle definitions.
Triangle 1 sets the global coordinate system:
- Node 1 at (0, 0)
- Node 2 at (L12, 0) (along positive X-axis)
- Node 3 via Law of Cosines
Subsequent triangles are attached to their shared edge.
"""
function reconstruct_geometry_from_triangles(num_nodes::Int, triangles::Vector{TriangleDef})
    coords = Dict{Int, Vector{Float64}}()
    triangle_history = Vector{Dict{String, Any}}()
    
    if isempty(triangles)
        return (coords, triangle_history)
    end
    
    # Process Triangle 1: Reference Triangle
    tri1 = triangles[1]
    n1, n2, n3 = tri1.nodes[1], tri1.nodes[2], tri1.nodes[3]
    L12, L23, L31 = tri1.lengths[1], tri1.lengths[2], tri1.lengths[3]
    
    valid1, msg1 = validate_triangle_lengths(L12, L23, L31)
    if !valid1
        error("Triangle 1 is invalid: $msg1")
    end
    
    # Place Node 1 and Node 2 along X-axis
    coords[n1] = [0.0, 0.0]
    coords[n2] = [L12, 0.0]
    
    # Calculate Node 3
    # c = L12, b = L31 (N1-N3), a = L23 (N2-N3)
    cos_theta1 = clamp((L12^2 + L31^2 - L23^2) / (2.0 * L12 * L31), -1.0, 1.0)
    theta1 = acos(cos_theta1)
    
    y_sign = (tri1.side == "right") ? -1.0 : 1.0
    x3 = L31 * cos_theta1
    y3 = y_sign * L31 * sin(theta1)
    coords[n3] = [x3, y3]
    
    push!(triangle_history, Dict{String, Any}(
        "triangle_id" => tri1.id,
        "nodes" => [n1, n2, n3],
        "shared_edge" => nothing,
        "new_node" => n3,
        "placement_side" => (y_sign > 0 ? "left" : "right")
    ))
    
    # Process subsequent triangles
    for k in 2:length(triangles)
        tri = triangles[k]
        t_nodes = tri.nodes
        t_lens = tri.lengths
        
        valid_k, msg_k = validate_triangle_lengths(t_lens[1], t_lens[2], t_lens[3])
        if !valid_k
            error("Triangle $(tri.id) is invalid: $msg_k")
        end
        
        known_nodes = [n for n in t_nodes if haskey(coords, n)]
        
        if length(known_nodes) < 2
            error("Triangle $(tri.id) ($(t_nodes)) does not share an edge with previously constructed geometry. Cannot place.")
        end
        
        if length(known_nodes) == 3
            # Closure triangle (all nodes already placed, adds bracing elements)
            push!(triangle_history, Dict{String, Any}(
                "triangle_id" => tri.id,
                "nodes" => t_nodes,
                "shared_edge" => "all_nodes_exist",
                "new_node" => nothing,
                "placement_side" => "internal"
            ))
            continue
        end
        
        # Exactly 2 known nodes, 1 new node
        A, B, C, c_user, b_len, a_len = find_shared_edge_info(t_nodes, t_lens, coords)
        
        pA = coords[A]
        pB = coords[B]
        
        # Calculate actual shared edge length from coordinates
        L_actual = hypot(pB[1] - pA[1], pB[2] - pA[2])
        c_len = L_actual # Use actual length for exact geometric closure
        
        # Validate triangle inequality with actual shared edge length
        val_sub, msg_sub = validate_triangle_lengths(a_len, b_len, c_len)
        if !val_sub
            error("Triangle $(tri.id) cannot be attached to shared edge $A-$B (length $(round(c_len, digits=4))): $msg_sub")
        end
        
        # Law of Cosines at Node A in local coordinate system
        cos_alpha = clamp((c_len^2 + b_len^2 - a_len^2) / (2.0 * c_len * b_len), -1.0, 1.0)
        alpha = acos(cos_alpha)
        
        # Edge orientation of A -> B in global coordinates
        theta_AB = atan(pB[2] - pA[2], pB[1] - pA[1])
        
        # Candidate positions: Left (+alpha) and Right (-alpha)
        x_loc = b_len * cos_alpha
        y_loc_pos = b_len * sin(alpha)
        
        pC_left = [
            pA[1] + x_loc * cos(theta_AB) - y_loc_pos * sin(theta_AB),
            pA[2] + x_loc * sin(theta_AB) + y_loc_pos * cos(theta_AB)
        ]
        
        pC_right = [
            pA[1] + x_loc * cos(theta_AB) - (-y_loc_pos) * sin(theta_AB),
            pA[2] + x_loc * sin(theta_AB) + (-y_loc_pos) * cos(theta_AB)
        ]
        
        # Decide placement side
        chosen_side = tri.side
        if chosen_side == "auto"
            # Automatic heuristic: choose side that does not overlap existing triangles,
            # or choose the candidate pointing outward from the existing truss centroid
            existing_centroid = [0.0, 0.0]
            for (_, pt) in coords
                existing_centroid .+= pt
            end
            existing_centroid ./= length(coords)
            
            d_left = hypot(pC_left[1] - existing_centroid[1], pC_left[2] - existing_centroid[2])
            d_right = hypot(pC_right[1] - existing_centroid[1], pC_right[2] - existing_centroid[2])
            
            # Prefer pointing outward away from centroid
            chosen_side = (d_left >= d_right) ? "left" : "right"
        end
        
        coords[C] = (chosen_side == "left") ? pC_left : pC_right
        
        push!(triangle_history, Dict{String, Any}(
            "triangle_id" => tri.id,
            "nodes" => t_nodes,
            "shared_edge" => [A, B],
            "new_node" => C,
            "placement_side" => chosen_side
        ))
    end
    
    return (coords, triangle_history)
end

end # module Triangles
