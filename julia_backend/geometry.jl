# julia_backend/geometry.jl
# Geometric calculation functions: Law of Cosines, coordinate transforms, angles

module Geometry

export calculate_triangle_angles, calculate_edge_orientation, transform_local_to_global,
       distance_2d, triangle_centroid, point_in_triangle, signed_triangle_area

"""
    calculate_triangle_angles(a::Float64, b::Float64, c::Float64) -> (alpha::Float64, beta::Float64, gamma::Float64)

Calculates the three interior angles of a triangle with side lengths (a, b, c) using the Law of Cosines.
- alpha is opposite to side a
- beta is opposite to side b
- gamma is opposite to side c
Angles are returned in radians.
"""
function calculate_triangle_angles(a::Float64, b::Float64, c::Float64)
    # Cosine law: cos(alpha) = (b^2 + c^2 - a^2) / (2*b*c)
    cos_a = clamp((b^2 + c^2 - a^2) / (2.0 * b * c), -1.0, 1.0)
    cos_b = clamp((a^2 + c^2 - b^2) / (2.0 * a * c), -1.0, 1.0)
    cos_c = clamp((a^2 + b^2 - c^2) / (2.0 * a * b), -1.0, 1.0)
    
    alpha = acos(cos_a)
    beta = acos(cos_b)
    gamma = acos(cos_c)
    
    return (alpha, beta, gamma)
end

"""
    calculate_edge_orientation(x1::Float64, y1::Float64, x2::Float64, y2::Float64) -> Float64

Returns the orientation angle theta (in radians) of the directed segment from (x1, y1) to (x2, y2).
Range: [-pi, pi].
"""
function calculate_edge_orientation(x1::Float64, y1::Float64, x2::Float64, y2::Float64)
    return atan(y2 - y1, x2 - x1)
end

"""
    transform_local_to_global(x_loc::Float64, y_loc::Float64, x_ref::Float64, y_ref::Float64, theta::Float64) -> (Float64, Float64)

Transforms a 2D point from a local coordinate system (with origin at (x_ref, y_ref) and rotated by theta)
into the global coordinate system:
    x_glob = x_ref + x_loc * cos(theta) - y_loc * sin(theta)
    y_glob = y_ref + x_loc * sin(theta) + y_loc * cos(theta)
"""
function transform_local_to_global(x_loc::Float64, y_loc::Float64, x_ref::Float64, y_ref::Float64, theta::Float64)
    c = cos(theta)
    s = sin(theta)
    x_glob = x_ref + x_loc * c - y_loc * s
    y_glob = y_ref + x_loc * s + y_loc * c
    return (x_glob, y_glob)
end

"""
    distance_2d(p1::Vector{Float64}, p2::Vector{Float64}) -> Float64
"""
function distance_2d(p1::Vector{Float64}, p2::Vector{Float64})
    return hypot(p2[1] - p1[1], p2[2] - p1[2])
end

"""
    triangle_centroid(p1::Vector{Float64}, p2::Vector{Float64}, p3::Vector{Float64}) -> Vector{Float64}
"""
function triangle_centroid(p1::Vector{Float64}, p2::Vector{Float64}, p3::Vector{Float64})
    return [(p1[1] + p2[1] + p3[1]) / 3.0, (p1[2] + p2[2] + p3[2]) / 3.0]
end

"""
    signed_triangle_area(p1::Vector{Float64}, p2::Vector{Float64}, p3::Vector{Float64}) -> Float64

Calculates 0.5 * cross_product(p2-p1, p3-p1). Positive if (p1, p2, p3) are counterclockwise (left turn).
"""
function signed_triangle_area(p1::Vector{Float64}, p2::Vector{Float64}, p3::Vector{Float64})
    return 0.5 * ((p2[1] - p1[1]) * (p3[2] - p1[2]) - (p2[2] - p1[2]) * (p3[1] - p1[1]))
end

"""
    point_in_triangle(pt::Vector{Float64}, p1::Vector{Float64}, p2::Vector{Float64}, p3::Vector{Float64}) -> Bool
"""
function point_in_triangle(pt::Vector{Float64}, p1::Vector{Float64}, p2::Vector{Float64}, p3::Vector{Float64})
    a1 = signed_triangle_area(pt, p1, p2)
    a2 = signed_triangle_area(pt, p2, p3)
    a3 = signed_triangle_area(pt, p3, p1)
    has_neg = (a1 < -1e-8) || (a2 < -1e-8) || (a3 < -1e-8)
    has_pos = (a1 > 1e-8) || (a2 > 1e-8) || (a3 > 1e-8)
    return !(has_neg && has_pos)
end

end # module Geometry
