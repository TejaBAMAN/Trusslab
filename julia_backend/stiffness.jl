# julia_backend/stiffness.jl
# Direct Stiffness element formulation for 2D truss members

module Stiffness

using LinearAlgebra

export calculate_element_stiffness, calculate_local_stiffness

"""
    calculate_element_stiffness(E::Float64, A::Float64, L::Float64, c::Float64, s::Float64) -> Matrix{Float64}

Computes the 4x4 element stiffness matrix in global coordinates:
    k_e = (A * E / L) * [
         c^2    c*s   -c^2   -c*s
         c*s    s^2   -c*s   -s^2
        -c^2   -c*s    c^2    c*s
        -c*s   -s^2    c*s    s^2
    ]
Degrees of freedom order: [u_ix, u_iy, u_jx, u_jy]
"""
function calculate_element_stiffness(E::Float64, A::Float64, L::Float64, c::Float64, s::Float64)
    k0 = (A * E) / L
    
    c2 = c * c
    s2 = s * s
    cs = c * s
    
    ke = zeros(Float64, 4, 4)
    
    ke[1, 1] =  c2
    ke[1, 2] =  cs
    ke[1, 3] = -c2
    ke[1, 4] = -cs
    
    ke[2, 1] =  cs
    ke[2, 2] =  s2
    ke[2, 3] = -cs
    ke[2, 4] = -s2
    
    ke[3, 1] = -c2
    ke[3, 2] = -cs
    ke[3, 3] =  c2
    ke[3, 4] =  cs
    
    ke[4, 1] = -cs
    ke[4, 2] = -s2
    ke[4, 3] =  cs
    ke[4, 4] =  s2
    
    return k0 .* ke
end

"""
    calculate_local_stiffness(E::Float64, A::Float64, L::Float64) -> Matrix{Float64}

Computes the 2x2 local element axial stiffness:
    k_local = (A * E / L) * [ 1 -1; -1 1 ]
"""
function calculate_local_stiffness(E::Float64, A::Float64, L::Float64)
    k0 = (A * E) / L
    return k0 .* [1.0 -1.0; -1.0 1.0]
end

end # module Stiffness
