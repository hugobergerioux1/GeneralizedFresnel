# Mesh, integration measures and mixed finite element spaces.

orderL = order          # Lagrange elements (Ey, and the multiplier of the SP formulation)
orderN = order - 1      # Nédélec elements (Ez, Ex)
degree = 2 * orderL     # quadrature degree

if mesh == "T"
    # triangles: Gmsh mesh, periodic in x, with the physical groups of the ports
    println("\nLoading mesh $meshfile ...")
    model = GmshDiscreteModel(meshfile)
    left_b, right_b = extrema(x -> x[1], get_node_coordinates(model))
    bot_b, top_b = extrema(x -> x[2], get_node_coordinates(model))
    port_tags_L = ["LportEdges","LportNodes"]
    port_tags_R = ["RportEdges","RportNodes"]
elseif mesh == "R"
    # rectangles: structured grid of 5 × 5 wavelengths of the incident wave, periodic in x
    println("\nBuilding the Cartesian grid $(partition[1]) × $(partition[2]) ...")
    λz = 2π / (k₀*real(sol.n[1,3]))
    λx = nx == 0 ? 2π / (√2*k₀) : 2π / (nx*k₀)
    left_b, right_b = -5λz/2, 5λz/2
    bot_b, top_b = -5λx/2, 5λx/2
    model = CartesianDiscreteModel((left_b, right_b, bot_b, top_b), partition; isperiodic=(false, true))
    iseven(partition[1]) || @warn "Odd number of cells along z: the interface z = 0 is not a grid line."

    # tags of the ports: the edges of the left and right sides
    topology = get_grid_topology(model)
    labels = get_face_labeling(model)
    for i in 1:num_faces(topology, 1)
        local edge_coords = get_face_coordinates(topology, 1)[i]
        local zmid = (edge_coords[1][1] + edge_coords[2][1])/2      # z of the middle of the edge
        if abs(zmid-left_b) < 1e-8
            model.face_labeling.d_to_dface_to_entity[2][i] = 3  # left wall
        elseif abs(zmid-right_b) < 1e-8
            model.face_labeling.d_to_dface_to_entity[2][i] = 4  # right wall
        end
    end
    # and their vertices
    vertex_coords = get_vertex_coordinates(topology)
    for i in 1:num_faces(topology, 0)
        if abs(vertex_coords[i][1]-left_b) < 1e-8
            model.face_labeling.d_to_dface_to_entity[1][i] = 3
        elseif abs(vertex_coords[i][1]-right_b) < 1e-8
            model.face_labeling.d_to_dface_to_entity[1][i] = 4
        end
    end
    add_tag_from_tags!(labels, "LportEdges", 3)
    add_tag_from_tags!(labels, "RportEdges", 4)
    port_tags_L = ["LportEdges"]
    port_tags_R = ["RportEdges"]
else
    error("mesh must be \"T\" (triangles) or \"R\" (rectangles)")
end

Ω = Triangulation(model)
Px = get_physical_coordinate(Ω)
dΩ = Gridap.Measure(Ω,degree)

# ports: left (z = left_b) and right (z = right_b) sides; top and bottom sides are periodic
Γₗ = BoundaryTriangulation(model, tags=port_tags_L)
dΓₗ = Gridap.Measure(Γₗ,degree)
Γᵣ = BoundaryTriangulation(model, tags=port_tags_R)
dΓᵣ = Gridap.Measure(Γᵣ,degree)

# domain and effective mesh size
Lz, Lx = right_b - left_b, top_b - bot_b
nelm = num_cells(Ω)
heff = mesh == "T" ? sqrt(4*Lz*Lx/(nelm*√3)) : sqrt(Lz*Lx/nelm)

λx = 2π / (k₀*nx)
if nx != 0 && abs(Lx/λx - round(Lx/λx)) > 1e-6
    @warn "The mesh width Lx = $Lx m is not a multiple of λx = $λx m: the periodic boundary " *
          "conditions are not compatible with the incident wave (check freq and nx of the case)."
end

# mesh with its physical groups, for visualization
mkpath(outdir)
writevtk(model, joinpath(outdir, "mesh_$(mesh)_$(round(Int, nelm/1000))k"))

# FE spaces: (Ez, Ex) in H(curl), Ey in H¹
V = FESpace(model, ReferenceFE(nedelec,Float64,orderN),conformity=:HCurl, vector_type=Vector{ComplexF64})
Q = FESpace(model, ReferenceFE(lagrangian,Float64,orderL),conformity=:H1, vector_type=Vector{ComplexF64})
U = TrialFESpace(V)
P = TrialFESpace(Q)
if formulation == "W"
    Y = MultiFieldFESpace([V, Q])           # test space
    X = MultiFieldFESpace([U, P])           # trial space
elseif formulation == "SP"
    T = FESpace(model, ReferenceFE(lagrangian,Float64,orderL),conformity=:H1, vector_type=Vector{ComplexF64},
        dirichlet_tags = [port_tags_L; port_tags_R])
    S = TrialFESpace(T,0)
    Y = MultiFieldFESpace([V, Q, T])        # test space
    X = MultiFieldFESpace([U, P, S])        # trial space
else
    error("formulation must be \"W\" (wave equation) or \"SP\" (saddle point)")
end
ndofs = num_free_dofs(X)

@printf("Mesh: %d %s, %d DOFs, h_eff = %.4g m, k₀ h_eff = %.4f, domain %.4g m × %.4g m\n",
        nelm, mesh == "T" ? "triangles" : "rectangles", ndofs, heff, k₀*heff, Lz, Lx)
println("Formulation $formulation, numerical collisionality Z = $Z")
