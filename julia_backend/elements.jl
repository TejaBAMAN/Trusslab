# julia_backend/elements.jl
# Automatic unique element generation, lengths, angles, direction cosines,
# and complete 9-shape cross-sectional property calculations & buckling configurations.

module Elements

using LinearAlgebra
using Printf

export Element, generate_unique_elements,
       SectionProperties, BucklingRestraint, MemberConfig,
       compute_solid_circle, compute_hollow_circle, compute_solid_rectangle, compute_solid_square,
       compute_rhs, compute_shs, compute_i_section, compute_tee_section,
       compute_channel_section, compute_angle_section, compute_custom_section,
       compute_section_properties

# =============================================================================
# Structural Element Definition
# =============================================================================
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

# =============================================================================
# Cross-Section Properties Data Structure
# =============================================================================
struct SectionProperties
    name::String
    shape::String
    params::Dict{String, Float64} # in SI units (m)
    A::Float64                    # Cross-sectional area (m^2)
    Ix::Float64                   # Centroidal moment of inertia about x-axis (m^4)
    Iy::Float64                   # Centroidal moment of inertia about y-axis (m^4)
    Ixy::Float64                  # Product of inertia about centroidal x-y axes (m^4)
    Iu::Float64                   # Principal second moment of area (u-axis) (m^4)
    Iv::Float64                   # Principal second moment of area (v-axis) (m^4)
    Imin::Float64                 # Minimum principal second moment of area (m^4)
    Imax::Float64                 # Maximum principal second moment of area (m^4)
    theta_p_rad::Float64          # Principal angle (rad)
    x_bar::Float64                # Centroid offset from reference origin (m)
    y_bar::Float64                # Centroid offset from reference origin (m)
end

# =============================================================================
# Buckling Restraints & Unbraced Lengths Data Structure
# =============================================================================
struct BucklingRestraint
    preset_u::String # "pinned_pinned", "fixed_fixed", "fixed_pinned", "cantilever", "custom"
    preset_v::String
    ku::Float64      # Effective length factor about u-axis (> 0)
    kv::Float64      # Effective length factor about v-axis (> 0)
    Lu::Float64      # Unbraced length about u-axis (> 0, m)
    Lv::Float64      # Unbraced length about v-axis (> 0, m)
    configured::Bool # True if verified by user, false if legacy/unconfigured
end

# Default constructor for unconfigured legacy members
BucklingRestraint(L::Float64) = BucklingRestraint("pinned_pinned", "pinned_pinned", 1.0, 1.0, L, L, false)

# =============================================================================
# Member Configuration (Combines Element, Section, Buckling, Material)
# =============================================================================
struct MemberConfig
    element_id::Int
    section::SectionProperties
    buckling::BucklingRestraint
    E::Float64       # Pa
    density::Float64 # kg/m^3
end

# =============================================================================
# Section Geometry Engine: 9 Built-in Shapes + Custom
# All geometric functions take SI dimensions in METRES.
# =============================================================================

"""
    compute_solid_circle(D::Float64; name="Solid Circle") -> SectionProperties
Solid circular rod of diameter D (m).
A = pi * D^2 / 4, Ix = Iy = pi * D^4 / 64.
"""
function compute_solid_circle(D::Float64; name::String="Solid Circle")
    D > 0 || error("Solid circle diameter D must be positive, got $D")
    A = pi * (D^2) / 4.0
    I = pi * (D^4) / 64.0
    params = Dict("D" => D)
    return SectionProperties(name, "solid_circle", params, A, I, I, 0.0, I, I, I, I, 0.0, D/2.0, D/2.0)
end

"""
    compute_hollow_circle(D::Float64, d_or_t::Float64; is_thickness::Bool=false, name="Circular Hollow Tube") -> SectionProperties
Circular hollow section (CHS / tube / pipe) with outer diameter D and inner diameter d (or wall t).
A = pi * (D^2 - d^2) / 4, Ix = Iy = pi * (D^4 - d^4) / 64.
"""
function compute_hollow_circle(D::Float64, d_or_t::Float64; is_thickness::Bool=false, name::String="Circular Hollow Tube")
    D > 0 || error("Circular tube outer diameter D must be positive, got $D")
    d = is_thickness ? (D - 2.0 * d_or_t) : d_or_t
    t = is_thickness ? d_or_t : (D - d) / 2.0
    
    t > 0 || error("Circular tube wall thickness must be positive, got $t")
    d > 0 || error("Circular tube inner diameter d must be positive, got $d (D=$D, t=$t)")
    D > d || error("Outer diameter D ($D) must exceed inner diameter d ($d)")

    A = pi * (D^2 - d^2) / 4.0
    I = pi * (D^4 - d^4) / 64.0
    params = Dict("D" => D, "d" => d, "t" => t)
    return SectionProperties(name, "hollow_circle", params, A, I, I, 0.0, I, I, I, I, 0.0, D/2.0, D/2.0)
end

"""
    compute_solid_rectangle(b::Float64, h::Float64; name="Solid Rectangle") -> SectionProperties
Solid rectangle of width b and depth h (m).
A = b * h, Ix = b * h^3 / 12, Iy = h * b^3 / 12.
"""
function compute_solid_rectangle(b::Float64, h::Float64; name::String="Solid Rectangle")
    b > 0 || error("Rectangle width b must be positive, got $b")
    h > 0 || error("Rectangle depth h must be positive, got $h")
    A = b * h
    Ix = (b * (h^3)) / 12.0
    Iy = (h * (b^3)) / 12.0
    Imin = min(Ix, Iy)
    Imax = max(Ix, Iy)
    params = Dict("b" => b, "h" => h)
    return SectionProperties(name, "solid_rectangle", params, A, Ix, Iy, 0.0, Ix, Iy, Imin, Imax, 0.0, b/2.0, h/2.0)
end

"""
    compute_solid_square(b::Float64; name="Solid Square") -> SectionProperties
Solid square with side b (m).
"""
function compute_solid_square(b::Float64; name::String="Solid Square")
    return compute_solid_rectangle(b, b; name=name)
end

"""
    compute_rhs(B::Float64, H::Float64, t::Float64; name="Rectangular Hollow Section (RHS)") -> SectionProperties
Rectangular hollow section (RHS) with outer dimensions B x H and uniform wall thickness t (m).
Ideal sharp-corner geometry:
b = B - 2t, h = H - 2t
A = B*H - b*h, Ix = (B*H^3 - b*h^3)/12, Iy = (H*B^3 - h*b^3)/12.
"""
function compute_rhs(B::Float64, H::Float64, t::Float64; name::String="Rectangular Hollow Section (RHS)")
    B > 0 || error("RHS outer width B must be positive, got $B")
    H > 0 || error("RHS outer depth H must be positive, got $H")
    t > 0 || error("RHS wall thickness t must be positive, got $t")
    2.0 * t < B || error("RHS wall thickness 2t ($((2*t))) must be less than outer width B ($B)")
    2.0 * t < H || error("RHS wall thickness 2t ($((2*t))) must be less than outer depth H ($H)")

    b = B - 2.0 * t
    h = H - 2.0 * t
    A = B * H - b * h
    Ix = (B * (H^3) - b * (h^3)) / 12.0
    Iy = (H * (B^3) - h * (b^3)) / 12.0
    Imin = min(Ix, Iy)
    Imax = max(Ix, Iy)
    params = Dict("B" => B, "H" => H, "t" => t, "b" => b, "h" => h)
    return SectionProperties(name, "rectangular_hollow", params, A, Ix, Iy, 0.0, Ix, Iy, Imin, Imax, 0.0, B/2.0, H/2.0)
end

"""
    compute_shs(B::Float64, t::Float64; name="Square Hollow Section (SHS)") -> SectionProperties
Square hollow section with outer dimension B and uniform wall thickness t.
"""
function compute_shs(B::Float64, t::Float64; name::String="Square Hollow Section (SHS)")
    return compute_rhs(B, B, t; name=name)
end

"""
    compute_i_section(B::Float64, H::Float64, tw::Float64, tf::Float64; name="Symmetric I-Section") -> SectionProperties
Symmetric I/H section, overall flange width B, depth H, web thickness tw, flange thickness tf.
hw = H - 2*tf
A = 2*B*tf + tw*hw
Ix = (B*H^3 - (B - tw)*hw^3) / 12
Iy = (2*tf*B^3 + hw*tw^3) / 12.
"""
function compute_i_section(B::Float64, H::Float64, tw::Float64, tf::Float64; name::String="Symmetric I-Section")
    B > 0 || error("I-section flange width B must be positive, got $B")
    H > 0 || error("I-section total depth H must be positive, got $H")
    tw > 0 || error("I-section web thickness tw must be positive, got $tw")
    tf > 0 || error("I-section flange thickness tf must be positive, got $tf")
    H > 2.0 * tf || error("I-section depth H ($H) must be greater than 2*tf ($((2*tf)))")
    B > tw || error("I-section flange width B ($B) must be greater than web thickness tw ($tw)")

    hw = H - 2.0 * tf
    A = 2.0 * B * tf + tw * hw
    Ix = (B * (H^3) - (B - tw) * (hw^3)) / 12.0
    Iy = (2.0 * tf * (B^3) + hw * (tw^3)) / 12.0
    Imin = min(Ix, Iy)
    Imax = max(Ix, Iy)
    params = Dict("B" => B, "H" => H, "tw" => tw, "tf" => tf, "hw" => hw)
    return SectionProperties(name, "i_section", params, A, Ix, Iy, 0.0, Ix, Iy, Imin, Imax, 0.0, B/2.0, H/2.0)
end

"""
    compute_tee_section(B::Float64, H::Float64, tw::Float64, tf::Float64; name="T-Section") -> SectionProperties
Built-up T-section with top flange B x tf and vertical web tw x (H - tf).
Uses parallel-axis theorem on non-overlapping rectangles.
"""
function compute_tee_section(B::Float64, H::Float64, tw::Float64, tf::Float64; name::String="T-Section")
    B > 0 || error("T-section flange width B must be positive, got $B")
    H > 0 || error("T-section total depth H must be positive, got $H")
    tw > 0 || error("T-section web thickness tw must be positive, got $tw")
    tf > 0 || error("T-section flange thickness tf must be positive, got $tf")
    H > tf || error("T-section depth H ($H) must exceed flange thickness tf ($tf)")
    B >= tw || error("T-section flange width B ($B) must be at least web thickness tw ($tw)")

    hw = H - tf
    # Rectangle 1: Top flange (B x tf)
    A1 = B * tf
    x1 = 0.0
    y1 = H - tf / 2.0
    Ix1 = B * (tf^3) / 12.0
    Iy1 = tf * (B^3) / 12.0

    # Rectangle 2: Web (tw x hw)
    A2 = tw * hw
    x2 = 0.0
    y2 = hw / 2.0
    Ix2 = tw * (hw^3) / 12.0
    Iy2 = hw * (tw^3) / 12.0

    A = A1 + A2
    x_bar = 0.0 # symmetric about y
    y_bar = (A1 * y1 + A2 * y2) / A

    Ix = Ix1 + A1 * ((y1 - y_bar)^2) + Ix2 + A2 * ((y2 - y_bar)^2)
    Iy = Iy1 + Iy2
    Ixy = 0.0 # symmetric about vertical axis
    Imin = min(Ix, Iy)
    Imax = max(Ix, Iy)

    params = Dict("B" => B, "H" => H, "tw" => tw, "tf" => tf, "hw" => hw)
    return SectionProperties(name, "tee_section", params, A, Ix, Iy, 0.0, Ix, Iy, Imin, Imax, 0.0, x_bar, y_bar)
end

"""
    compute_channel_section(B::Float64, H::Float64, tw::Float64, tf::Float64; name="Channel Section") -> SectionProperties
Built-up C-channel with web tw x (H - 2*tf) and top/bottom flanges B x tf.
Parallel-axis theorem on non-overlapping rectangles.
"""
function compute_channel_section(B::Float64, H::Float64, tw::Float64, tf::Float64; name::String="Channel Section")
    B > 0 || error("Channel flange width B must be positive, got $B")
    H > 0 || error("Channel depth H must be positive, got $H")
    tw > 0 || error("Channel web thickness tw must be positive, got $tw")
    tf > 0 || error("Channel flange thickness tf must be positive, got $tf")
    H > 2.0 * tf || error("Channel depth H ($H) must exceed 2*tf ($((2*tf)))")
    B > tw || error("Channel flange width B ($B) must exceed web thickness tw ($tw)")

    hw = H - 2.0 * tf
    Aw = tw * hw
    xw = tw / 2.0
    yw = H / 2.0
    Ixw = tw * (hw^3) / 12.0
    Iyw = hw * (tw^3) / 12.0

    Af = B * tf
    xft = B / 2.0
    yft = H - tf / 2.0
    xfb = B / 2.0
    yfb = tf / 2.0
    Ixf = B * (tf^3) / 12.0
    Iyf = tf * (B^3) / 12.0

    A = Aw + 2.0 * Af
    x_bar = (Aw * xw + 2.0 * Af * xft) / A
    y_bar = H / 2.0 # symmetric horizontally

    Ix = Ixw + Aw * ((yw - y_bar)^2) + 2.0 * (Ixf + Af * ((yft - y_bar)^2))
    Iy = Iyw + Aw * ((xw - x_bar)^2) + 2.0 * (Iyf + Af * ((xft - x_bar)^2))
    Ixy = 0.0 # horizontal symmetry implies zero product of inertia
    Imin = min(Ix, Iy)
    Imax = max(Ix, Iy)

    params = Dict("B" => B, "H" => H, "tw" => tw, "tf" => tf, "hw" => hw)
    return SectionProperties(name, "channel_section", params, A, Ix, Iy, 0.0, Ix, Iy, Imin, Imax, 0.0, x_bar, y_bar)
end

"""
    compute_angle_section(B::Float64, H::Float64, th::Float64, tv::Float64; name="Angle / L-Section") -> SectionProperties
Built-up equal or unequal L-angle section.
Vertical leg: tv x (H - th), Horizontal leg: B x th.
Computes area, centroid, centroidal Ix, Iy, Ixy and exact principal second moments (Imin, Imax).
Note: For unsymmetric angles, Ixy is nonzero, so Imin < min(Ix, Iy).
"""
function compute_angle_section(B::Float64, H::Float64, th::Float64, tv::Float64; name::String="Angle Section")
    B > 0 || error("Angle width B must be positive, got $B")
    H > 0 || error("Angle height H must be positive, got $H")
    th > 0 || error("Horizontal leg thickness th must be positive, got $th")
    tv > 0 || error("Vertical leg thickness tv must be positive, got $tv")
    B > tv || error("Angle horizontal leg B ($B) must exceed vertical leg thickness tv ($tv)")
    H > th || error("Angle vertical leg H ($H) must exceed horizontal leg thickness th ($th)")

    # Component 1: vertical leg above horizontal leg: tv x (H - th)
    h1 = H - th
    A1 = tv * h1
    x1 = tv / 2.0
    y1 = th + h1 / 2.0
    Ix1 = tv * (h1^3) / 12.0
    Iy1 = h1 * (tv^3) / 12.0

    # Component 2: horizontal leg: B x th
    A2 = B * th
    x2 = B / 2.0
    y2 = th / 2.0
    Ix2 = B * (th^3) / 12.0
    Iy2 = th * (B^3) / 12.0

    A = A1 + A2
    x_bar = (A1 * x1 + A2 * x2) / A
    y_bar = (A1 * y1 + A2 * y2) / A

    Ix = Ix1 + A1 * ((y1 - y_bar)^2) + Ix2 + A2 * ((y2 - y_bar)^2)
    Iy = Iy1 + A1 * ((x1 - x_bar)^2) + Iy2 + A2 * ((x2 - x_bar)^2)
    Ixy = A1 * (x1 - x_bar) * (y1 - y_bar) + A2 * (x2 - x_bar) * (y2 - y_bar)

    # Principal second moments of area
    R = sqrt(((Ix - Iy) / 2.0)^2 + Ixy^2)
    Imin = (Ix + Iy) / 2.0 - R
    Imax = (Ix + Iy) / 2.0 + R
    theta_p = 0.5 * atan(-2.0 * Ixy, Ix - Iy)

    params = Dict("B" => B, "H" => H, "th" => th, "tv" => tv)
    return SectionProperties(name, "angle_section", params, A, Ix, Iy, Ixy, Imin, Imax, Imin, Imax, theta_p, x_bar, y_bar)
end

"""
    compute_custom_section(; A, Iu, Iv, Ix=0.0, Iy=0.0, Ixy=0.0, name="Custom") -> SectionProperties
Direct entry of custom / tabulated section properties.
Supports either (A, Iu, Iv) or (A, Ix, Iy, Ixy).
Validates physically admissible positive inertia tensor (Ix*Iy - Ixy^2 > 0).
"""
function compute_custom_section(;
    A::Float64,
    Iu::Float64=0.0,
    Iv::Float64=0.0,
    Ix::Float64=0.0,
    Iy::Float64=0.0,
    Ixy::Float64=0.0,
    name::String="Custom Section"
)
    A > 0 || error("Custom section area A must be positive, got $A")
    
    if Iu > 0 && Iv > 0 && Ix == 0.0 && Iy == 0.0
        # Direct principal second moments mode
        Imin = min(Iu, Iv)
        Imax = max(Iu, Iv)
        params = Dict("A" => A, "Iu" => Iu, "Iv" => Iv)
        return SectionProperties(name, "custom", params, A, Iu, Iv, 0.0, Iu, Iv, Imin, Imax, 0.0, 0.0, 0.0)
    else
        # Centroidal Ix, Iy, Ixy mode
        Ix > 0 || error("Custom section Ix must be positive, got $Ix")
        Iy > 0 || error("Custom section Iy must be positive, got $Iy")
        det = Ix * Iy - Ixy^2
        det > 0 || error("Custom section inertia tensor must be positive-definite (Ix*Iy - Ixy^2 > 0), got det=$det")

        R = sqrt(((Ix - Iy) / 2.0)^2 + Ixy^2)
        Imin = (Ix + Iy) / 2.0 - R
        Imax = (Ix + Iy) / 2.0 + R
        Imin > 0 || error("Computed minimum principal inertia must be positive, got $Imin")
        theta_p = 0.5 * atan(-2.0 * Ixy, Ix - Iy)
        params = Dict("A" => A, "Ix" => Ix, "Iy" => Iy, "Ixy" => Ixy)
        return SectionProperties(name, "custom", params, A, Ix, Iy, Ixy, Imin, Imax, Imin, Imax, theta_p, 0.0, 0.0)
    end
end

"""
    compute_section_properties(shape::String, params::Dict{String, Float64}; name::String="") -> SectionProperties
Dispatcher for all supported cross-section shapes.
"""
function compute_section_properties(shape::String, params::Dict{String, Float64}; name::String="")
    s = lowercase(strip(shape))
    default_name = isempty(name) ? uppercase(s) : name

    if s == "solid_circle" || s == "circle" || s == "rod"
        D = params["D"]
        return compute_solid_circle(D; name=default_name)
    elseif s == "hollow_circle" || s == "pipe" || s == "tube" || s == "chs"
        D = params["D"]
        if haskey(params, "t")
            return compute_hollow_circle(D, params["t"]; is_thickness=true, name=default_name)
        else
            return compute_hollow_circle(D, params["d"]; is_thickness=false, name=default_name)
        end
    elseif s == "solid_rectangle" || s == "rectangle"
        return compute_solid_rectangle(params["b"], params["h"]; name=default_name)
    elseif s == "solid_square" || s == "square"
        return compute_solid_square(params["b"]; name=default_name)
    elseif s == "rectangular_hollow" || s == "rhs"
        return compute_rhs(params["B"], params["H"], params["t"]; name=default_name)
    elseif s == "square_hollow" || s == "shs"
        return compute_shs(params["B"], params["t"]; name=default_name)
    elseif s == "i_section" || s == "i_beam" || s == "h_section" || s == "w_shape"
        return compute_i_section(params["B"], params["H"], params["tw"], params["tf"]; name=default_name)
    elseif s == "tee_section" || s == "t_section"
        return compute_tee_section(params["B"], params["H"], params["tw"], params["tf"]; name=default_name)
    elseif s == "channel_section" || s == "channel" || s == "c_section"
        return compute_channel_section(params["B"], params["H"], params["tw"], params["tf"]; name=default_name)
    elseif s == "angle_section" || s == "angle" || s == "l_section"
        th = get(params, "th", get(params, "t", 0.0))
        tv = get(params, "tv", get(params, "t", 0.0))
        return compute_angle_section(params["B"], params["H"], th, tv; name=default_name)
    elseif s == "custom"
        A = params["A"]
        Iu = get(params, "Iu", 0.0)
        Iv = get(params, "Iv", 0.0)
        Ix = get(params, "Ix", 0.0)
        Iy = get(params, "Iy", 0.0)
        Ixy = get(params, "Ixy", 0.0)
        return compute_custom_section(; A=A, Iu=Iu, Iv=Iv, Ix=Ix, Iy=Iy, Ixy=Ixy, name=default_name)
    else
        error("Unsupported section shape: '$shape'")
    end
end

# =============================================================================
# Automatic Unique Element Extraction from Triangles
# =============================================================================
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
