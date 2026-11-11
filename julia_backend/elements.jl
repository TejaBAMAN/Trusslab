# julia_backend/elements.jl
# Automatic unique element generation, lengths, angles, and direction cosines

module Elements

export Element, generate_unique_elements

struct Element
    id::Int
    node_i::Int
    node_j::Int
    length::Float64
    c::Float64        # cos(theta) = dx / L
    s::Float64        # sin(theta) = dy / L
    angle_rad::Float64
    angle_deg::Float64
end

"""
    generate_unique_elements(triangles::Vector{<:Any}, coords::Dict{Int, Vector{Float64}}) -> Vector{Element}

Extracts unique members from triangle definitions, ensuring no duplicate elements.
Calculates length, direction cosines, and angle for every element.
"""
function generate_unique_elements(triangles::Vector{<:Any}, coords::Dict{Int, Vector{Float64}})
    seen_edges = Set{Tuple{Int, Int}}()
    unique_pairs = Vector{Tuple{Int, Int}}()
    
    for tri in triangles
        nodes = tri.nodes
        n1, n2, n3 = nodes[1], nodes[2], nodes[3]
        
        edges = [
            (min(n1, n2), max(n1, n2)),
            (min(n2, n3), max(n2, n3)),
            (min(n3, n1), max(n3, n1))
        ]
        
        for e in edges
            if !(e in seen_edges)
                push!(seen_edges, e)
                push!(unique_pairs, e)
            end
        end
    end
    
    # Sort pairs for consistent canonical ordering
    sort!(unique_pairs)
    
    elements = Element[]
    for (idx, (ni, nj)) in enumerate(unique_pairs)
        if !haskey(coords, ni) || !haskey(coords, nj)
            error("Coordinates missing for node $ni or $nj in element $idx.")
        end
        
        pi = coords[ni]
        pj = coords[nj]
        
        dx = pj[1] - pi[1]
        dy = pj[2] - pi[2]
        L = hypot(dx, dy)
        
        if L < 1e-9
            error("Element $idx between node $ni and $nj has zero length.")
        end
        
        c = dx / L
        s = dy / L
        ang_rad = atan(dy, dx)
        ang_deg = rad2deg(ang_rad)
        
        push!(elements, Element(idx, ni, nj, L, c, s, ang_rad, ang_deg))
    end
    
    return elements
end

end # module Elements
