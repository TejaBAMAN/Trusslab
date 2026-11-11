# julia_backend/validation.jl
# Module for geometric and structural validation checks

module Validation

export validate_triangle_lengths, validate_connectivity, check_zero_length, check_duplicate_nodes

"""
    validate_triangle_lengths(a::Float64, b::Float64, c::Float64; tol::Float64=1e-7) -> (Bool, String)

Checks whether three side lengths can form a valid, non-degenerate triangle
using the Triangle Inequality Theorem:
    a + b > c
    b + c > a
    c + a > b
"""
function validate_triangle_lengths(a::Float64, b::Float64, c::Float64; tol::Float64=1e-7)
    if a <= tol || b <= tol || c <= tol
        return (false, "Side lengths must be strictly positive. Received ($a, $b, $c).")
    end
    
    if (a + b <= c + tol) || (b + c <= a + tol) || (c + a <= b + tol)
        return (false, "The specified side lengths ($a, $b, $c) cannot form a valid triangle because they violate the Triangle Inequality.")
    end
    
    return (true, "Valid triangle.")
end

"""
    validate_connectivity(num_nodes::Int, triangles::Vector{Any}) -> (Bool, String)

Validates that all node IDs referenced in triangles are within 1:num_nodes.
"""
function validate_connectivity(num_nodes::Int, triangles::Vector{<:Any})
    for (idx, tri) in enumerate(triangles)
        nodes = get(tri, "nodes", Int[])
        if length(nodes) != 3
            return (false, "Triangle $idx must specify exactly 3 node IDs.")
        end
        if length(unique(nodes)) != 3
            return (false, "Triangle $idx contains duplicate node IDs: $(nodes).")
        end
        for n in nodes
            if n < 1 || n > num_nodes
                return (false, "Triangle $idx references invalid node ID $n. Expected between 1 and $num_nodes.")
            end
        end
    end
    return (true, "Connectivity valid.")
end

"""
    check_zero_length(coords_i::Vector{Float64}, coords_j::Vector{Float64}; tol::Float64=1e-9) -> Bool
"""
function check_zero_length(coords_i::Vector{Float64}, coords_j::Vector{Float64}; tol::Float64=1e-9)
    dx = coords_j[1] - coords_i[1]
    dy = coords_j[2] - coords_i[2]
    return sqrt(dx*dx + dy*dy) < tol
end

"""
    check_duplicate_nodes(coords::Dict{Int, Vector{Float64}}; tol::Float64=1e-6) -> Vector{Tuple{Int, Int}}
"""
function check_duplicate_nodes(coords::Dict{Int, Vector{Float64}}; tol::Float64=1e-6)
    dups = Tuple{Int, Int}[]
    keys_arr = sort(collect(keys(coords)))
    for i in 1:length(keys_arr)
        for j in (i+1):length(keys_arr)
            n1 = keys_arr[i]
            n2 = keys_arr[j]
            p1 = coords[n1]
            p2 = coords[n2]
            if hypot(p1[1]-p2[1], p1[2]-p2[2]) < tol
                push!(dups, (n1, n2))
            end
        end
    end
    return dups
end

end # module Validation
