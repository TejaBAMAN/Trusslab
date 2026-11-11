# julia_backend/tests/runtests.jl
# Comprehensive automated test suite for Julia Truss Analysis Engine

using Test
using LinearAlgebra

# Include parent modules
parent_dir = dirname(@__DIR__)
include(joinpath(parent_dir, "validation.jl"))
include(joinpath(parent_dir, "geometry.jl"))
include(joinpath(parent_dir, "triangles.jl"))
include(joinpath(parent_dir, "elements.jl"))
include(joinpath(parent_dir, "stiffness.jl"))
include(joinpath(parent_dir, "boundary_conditions.jl"))
include(joinpath(parent_dir, "assembly.jl"))
include(joinpath(parent_dir, "stability.jl"))
include(joinpath(parent_dir, "results.jl"))
include(joinpath(parent_dir, "solver.jl"))
include(joinpath(parent_dir, "communication.jl"))

using .Validation
using .Geometry
using .Triangles
using .Elements
using .Stiffness
using .BoundaryConditions
using .Assembly
using .Stability
using .Results
using .Solver
using .Communication

println("==================================================")
println("RUNNING JULIA TRUSS ENGINE TEST SUITE (12 TESTS)")
println("==================================================")

@testset "Julia 2D Truss Engine Test Suite" begin

    # ----------------------------------------------------
    # Test 1: Triangle validity
    # ----------------------------------------------------
    @testset "Test 1: Triangle Validity" begin
        # Valid 3-4-5 triangle
        v1, msg1 = validate_triangle_lengths(3.0, 4.0, 5.0)
        @test v1 == true
        
        # Valid equilateral triangle
        v2, msg2 = validate_triangle_lengths(5.0, 5.0, 5.0)
        @test v2 == true
        
        # Invalid triangle: sum of two sides equals third (collinear)
        v3, msg3 = validate_triangle_lengths(2.0, 3.0, 5.0)
        @test v3 == false
        
        # Invalid triangle: one side greater than sum of others
        v4, msg4 = validate_triangle_lengths(1.0, 2.0, 10.0)
        @test v4 == false
        
        # Invalid: negative or zero length
        v5, msg5 = validate_triangle_lengths(-3.0, 4.0, 5.0)
        @test v5 == false
        v6, msg6 = validate_triangle_lengths(0.0, 4.0, 5.0)
        @test v6 == false
    end

    # ----------------------------------------------------
    # Test 2: Law of Cosines
    # ----------------------------------------------------
    @testset "Test 2: Law of Cosines" begin
        # 3-4-5 right triangle:
        # a = 3, b = 4, c = 5
        # alpha opposite a (3): asin(3/5) = 0.6435 rad (36.87 deg)
        # beta opposite b (4): asin(4/5) = 0.9273 rad (53.13 deg)
        # gamma opposite c (5): 90 deg = pi/2
        alpha, beta, gamma = calculate_triangle_angles(3.0, 4.0, 5.0)
        @test isapprox(rad2deg(gamma), 90.0; atol=1e-5)
        @test isapprox(rad2deg(alpha), 36.8698976; atol=1e-5)
        @test isapprox(rad2deg(beta), 53.1301024; atol=1e-5)
        @test isapprox(alpha + beta + gamma, pi; atol=1e-7)
        
        # Equilateral triangle: 60 deg each
        a_eq, b_eq, c_eq = calculate_triangle_angles(10.0, 10.0, 10.0)
        @test isapprox(rad2deg(a_eq), 60.0; atol=1e-5)
        @test isapprox(rad2deg(b_eq), 60.0; atol=1e-5)
        @test isapprox(rad2deg(c_eq), 60.0; atol=1e-5)
    end

    # ----------------------------------------------------
    # Test 3: Reference Triangle
    # ----------------------------------------------------
    @testset "Test 3: Reference Triangle" begin
        # Triangle 1 with nodes 1, 2, 3 and lengths L12=5.0, L23=3.0, L31=4.0
        # Reference edge 1-2 must be at 0 degrees along positive X-axis
        tri1 = TriangleDef(1, [1, 2, 3], [5.0, 3.0, 4.0], "left")
        coords, hist = reconstruct_geometry_from_triangles(3, [tri1])
        
        @test isapprox(coords[1], [0.0, 0.0]; atol=1e-7)
        @test isapprox(coords[2], [5.0, 0.0]; atol=1e-7)
        
        # Orientation of edge 1-2
        theta_12 = calculate_edge_orientation(coords[1][1], coords[1][2], coords[2][1], coords[2][2])
        @test isapprox(theta_12, 0.0; atol=1e-7)
    end

    # ----------------------------------------------------
    # Test 4: Node Coordinate Generation
    # ----------------------------------------------------
    @testset "Test 4: Node Coordinate Generation" begin
        # Triangle 1: L12=5.0, L23=3.0, L31=4.0
        # In right triangle 3-4-5 with hypotenuse 5 between nodes 1 and 2:
        # Node 1 at (0,0), Node 2 at (5,0)
        # Node 3 is at distance 4 from node 1 and 3 from node 2:
        # cos(theta1) = (5^2 + 4^2 - 3^2)/(2*5*4) = (25+16-9)/40 = 32/40 = 0.8
        # sin(theta1) = 0.6
        # x3 = 4 * 0.8 = 3.2
        # y3 = 4 * 0.6 = 2.4
        tri1 = TriangleDef(1, [1, 2, 3], [5.0, 3.0, 4.0], "left")
        coords, _ = reconstruct_geometry_from_triangles(3, [tri1])
        
        @test isapprox(coords[3][1], 3.2; atol=1e-5)
        @test isapprox(coords[3][2], 2.4; atol=1e-5)
        @test isapprox(hypot(coords[3][1], coords[3][2]), 4.0; atol=1e-5)
        @test isapprox(hypot(coords[3][1]-5.0, coords[3][2]), 3.0; atol=1e-5)
    end

    # ----------------------------------------------------
    # Test 5: Shared-Edge Detection
    # ----------------------------------------------------
    @testset "Test 5: Shared-Edge Detection" begin
        # Triangle 1: 1, 2, 3 (lengths 5, 3, 4)
        # Triangle 2: 2, 3, 4 (shared edge 2-3 with length 3, and 3-4=4, 4-2=5)
        tri1 = TriangleDef(1, [1, 2, 3], [5.0, 3.0, 4.0], "left")
        tri2 = TriangleDef(2, [2, 3, 4], [3.0, 4.0, 5.0], "right")
        
        coords, _ = reconstruct_geometry_from_triangles(4, [tri1, tri2])
        @test length(coords) == 4
        
        elements = generate_unique_elements([tri1, tri2], coords)
        # Should have exactly 5 elements (1-2, 2-3, 1-3, 3-4, 2-4), NOT 6 duplicate elements!
        @test length(elements) == 5
        
        # Verify edge 2-3 appears only once
        count_23 = count(e -> (e.node_i == 2 && e.node_j == 3), elements)
        @test count_23 == 1
    end

    # ----------------------------------------------------
    # Test 6: Tilted Triangle Construction
    # ----------------------------------------------------
    @testset "Test 6: Tilted Triangle Construction" begin
        # Edge 2-3 is tilted:
        # Node 2 is at (5.0, 0.0), Node 3 is at (3.2, 2.4)
        # Angle of 2-3 is atan2(2.4, -1.8) = ~126.87 deg (tilted!)
        # Triangle 2 is attached to tilted edge 2-3
        tri1 = TriangleDef(1, [1, 2, 3], [5.0, 3.0, 4.0], "left")
        tri2 = TriangleDef(2, [2, 3, 4], [3.0, 4.0, 5.0], "left")
        
        coords, hist = reconstruct_geometry_from_triangles(4, [tri1, tri2])
        p4 = coords[4]
        
        # Distance from node 3 to 4 must be 4.0
        d34 = hypot(p4[1] - coords[3][1], p4[2] - coords[3][2])
        @test isapprox(d34, 4.0; atol=1e-5)
        
        # Distance from node 2 to 4 must be 5.0
        d24 = hypot(p4[1] - coords[2][1], p4[2] - coords[2][2])
        @test isapprox(d24, 5.0; atol=1e-5)
    end

    # ----------------------------------------------------
    # Test 7: Element Length
    # ----------------------------------------------------
    @testset "Test 7: Element Length" begin
        tri1 = TriangleDef(1, [1, 2, 3], [6.0, 8.0, 10.0], "left")
        coords, _ = reconstruct_geometry_from_triangles(3, [tri1])
        elements = generate_unique_elements([tri1], coords)
        
        lengths = sort([e.length for e in elements])
        @test isapprox(lengths[1], 6.0; atol=1e-5)
        @test isapprox(lengths[2], 8.0; atol=1e-5)
        @test isapprox(lengths[3], 10.0; atol=1e-5)
    end

    # ----------------------------------------------------
    # Test 8: Direction Cosines
    # ----------------------------------------------------
    @testset "Test 8: Direction Cosines" begin
        tri1 = TriangleDef(1, [1, 2, 3], [5.0, 4.0, 3.0], "left")
        coords, _ = reconstruct_geometry_from_triangles(3, [tri1])
        elements = generate_unique_elements([tri1], coords)
        
        for e in elements
            @test isapprox(e.c^2 + e.s^2, 1.0; atol=1e-7)
        end
    end

    # ----------------------------------------------------
    # Test 9: Stiffness Matrix
    # ----------------------------------------------------
    @testset "Test 9: Stiffness Matrix" begin
        E = 200e9 # 200 GPa
        A = 0.01  # 0.01 m^2
        L = 4.0   # 4 m
        # Horizontal member (c=1, s=0)
        ke_horiz = calculate_element_stiffness(E, A, L, 1.0, 0.0)
        k0 = (A * E) / L # 5e8 N/m
        @test isapprox(ke_horiz[1, 1], k0; atol=1e-3)
        @test isapprox(ke_horiz[1, 3], -k0; atol=1e-3)
        @test isapprox(ke_horiz[2, 2], 0.0; atol=1e-3)
        @test isapprox(ke_horiz[4, 4], 0.0; atol=1e-3)
        
        # 45-degree member (c = 1/sqrt(2), s = 1/sqrt(2))
        c45 = 1.0 / sqrt(2.0)
        s45 = 1.0 / sqrt(2.0)
        ke_45 = calculate_element_stiffness(E, A, L, c45, s45)
        @test isapprox(ke_45[1, 1], 0.5 * k0; atol=1e-3)
        @test isapprox(ke_45[1, 2], 0.5 * k0; atol=1e-3)
        @test isapprox(ke_45[2, 2], 0.5 * k0; atol=1e-3)
        @test isapprox(ke_45[3, 3], 0.5 * k0; atol=1e-3)
        
        # Verify symmetry
        @test isapprox(ke_45, ke_45'; atol=1e-7)
    end

    # ----------------------------------------------------
    # Test 10: Simple Truss FEM Solution
    # ----------------------------------------------------
    @testset "Test 10: Simple Truss FEM Solution" begin
        # Classic symmetric 3-node triangular truss:
        # Base: Node 1 at (0,0), Node 2 at (4,0) -> L12 = 4
        # Node 3 at (2, 2*sqrt(3)): equilateral-like isosceles with L13 = 4, L23 = 4
        # Applied load at Node 3: Fy = -10,000 N
        # Pinned at Node 1 (fixed X, Y), Roller at Node 2 (fixed Y)
        # By symmetry:
        # Rx1 = 0
        # Ry1 = 5,000 N
        # Ry2 = 5,000 N
        # Node 3 angle with horizontal: 60 deg
        # Member 1-3 axial force: Compressive = 5000 / sin(60) = 5773.5 N (Compression)
        # Member 2-3 axial force: Compressive = 5000 / sin(60) = 5773.5 N (Compression)
        # Member 1-2 axial force: Tensile = 5000 / tan(60) = 2886.75 N (Tension)
        
        tri = TriangleDef(1, [1, 2, 3], [4.0, 4.0, 4.0], "left")
        coords, _ = reconstruct_geometry_from_triangles(3, [tri])
        elements = generate_unique_elements([tri], coords)
        
        E = 200e9
        A = 0.001 # 10 cm^2
        supports = [
            SupportDef(1, "pinned", true, true),
            SupportDef(2, "roller_x", false, true)
        ]
        loads = [LoadDef(3, 0.0, -10000.0)]
        
        res = solve_truss_fem(3, elements, coords, E, A, supports, loads)
        @test res.success == true
        
        # Check Reactions
        r_dict = Dict(r.node => r for r in res.reactions)
        @test isapprox(r_dict[1].rx, 0.0; atol=1e-3)
        @test isapprox(r_dict[1].ry, 5000.0; atol=1e-3)
        @test isapprox(r_dict[2].ry, 5000.0; atol=1e-3)
        
        # Check Element forces
        # Member 1-2 (horizontal base):
        el_12 = findfirst(e -> (e.node_i == 1 && e.node_j == 2), res.elements)
        @test isapprox(res.elements[el_12].axial_force, 2886.75; atol=1.0)
        @test res.elements[el_12].state == "Tension"
        
        # Member 1-3:
        el_13 = findfirst(e -> (e.node_i == 1 && e.node_j == 3), res.elements)
        @test isapprox(res.elements[el_13].axial_force, -5773.50; atol=1.0)
        @test res.elements[el_13].state == "Compression"
        
        # Member 2-3:
        el_23 = findfirst(e -> (e.node_i == 2 && e.node_j == 3), res.elements)
        @test isapprox(res.elements[el_23].axial_force, -5773.50; atol=1.0)
        @test res.elements[el_23].state == "Compression"
    end

    # ----------------------------------------------------
    # Test 11: Equilibrium Check
    # ----------------------------------------------------
    @testset "Test 11: Reactions and Equilibrium" begin
        # 4-node, 2-triangle truss with arbitrary loads: Fx=12000, Fy=-25000
        tri1 = TriangleDef(1, [1, 2, 3], [4.0, 3.0, 5.0], "left")
        tri2 = TriangleDef(2, [2, 3, 4], [3.0, 4.0, 5.0], "right")
        coords, _ = reconstruct_geometry_from_triangles(4, [tri1, tri2])
        elements = generate_unique_elements([tri1, tri2], coords)
        
        E = 210e9
        A = 0.002
        supports = [
            SupportDef(1, "pinned", true, true),
            SupportDef(4, "roller_x", false, true)
        ]
        loads = [
            LoadDef(2, 5000.0, -8000.0),
            LoadDef(3, 7000.0, -17000.0)
        ]
        
        res = solve_truss_fem(4, elements, coords, E, A, supports, loads)
        @test res.success == true
        
        total_fx = sum(ld.fx for ld in loads)
        total_fy = sum(ld.fy for ld in loads)
        sum_rx = sum(r.rx for r in res.reactions)
        sum_ry = sum(r.ry for r in res.reactions)
        
        @test isapprox(total_fx + sum_rx, 0.0; atol=1e-4)
        @test isapprox(total_fy + sum_ry, 0.0; atol=1e-4)
        @test res.summary.equilibrium_satisfied == true
    end

    # ----------------------------------------------------
    # Test 12: Unstable Truss Detection
    # ----------------------------------------------------
    @testset "Test 12: Unstable Truss Detection" begin
        # Triangle with only 2 support restraints (needs at least 3)
        tri1 = TriangleDef(1, [1, 2, 3], [5.0, 4.0, 3.0], "left")
        coords, _ = reconstruct_geometry_from_triangles(3, [tri1])
        elements = generate_unique_elements([tri1], coords)
        
        E = 200e9
        A = 0.001
        # Only 2 rollers in Y (free to translate horizontally and rotate!)
        insufficient_supports = [
            SupportDef(1, "roller_x", false, true),
            SupportDef(2, "roller_x", false, true)
        ]
        loads = [LoadDef(3, 0.0, -1000.0)]
        
        res_unstable = solve_truss_fem(3, elements, coords, E, A, insufficient_supports, loads)
        # Must detect instability and fail gracefully without crashing
        @test res_unstable.success == false
        @test contains(res_unstable.error_message, "Insufficient supports") || contains(res_unstable.error_message, "Rigid body")
    end

end

println("==================================================")
println("ALL 12 TESTS PASSED SUCCESSFULLY!")
println("==================================================")
