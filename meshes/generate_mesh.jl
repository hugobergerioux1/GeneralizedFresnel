# Meshes of the generalized Fresnel problem
#
#   julia --project=.. generate_mesh.jl
#
# regenerates the two reference meshes of the paper (incidence 45° from vacuum, 5.8 GHz):
#   box2D_6k.msh  (resol = 50,  6k triangles)
#   box2D_15k.msh (resol = 81, 15k triangles)
#
# change the value of resol for finer or coarser meshes

using GridapGmsh: gmsh         # Gmsh API shipped with GridapGmsh

"""
Rectangle z ∈ [−Lz/2, Lz/2] × x ∈ [−Lx/2, Lx/2], split by the interface z = 0, with
Lz = nwav λz and Lx = nwav λx, where λx = λ₀/nx and λz = λ₀/√(1 − nx²) are the wavelengths
along x and z of the incident wave (coming from vacuum). Triangles of characteristic size
lc = Lz/resol (Frontal-Delaunay). The top and bottom sides are periodic: since Lx is a
multiple of λx, the mesh is only valid for this frequency and this tangential index nx.

Physical groups: "LportEdges", "LportNodes" (left port z = −Lz/2), "RportEdges",
"RportNodes" (right port z = Lz/2), "Interface", "Medium 1" (z < 0), "Medium 2" (z > 0),
"top", "bottom".
"""
function generate_mesh(file; freq = 5.8e9, nx = sind(45), nwav = 5, resol = 81)
    ω = 2π * freq
    c = 299792458.0
    λx = nx == 0 ? 2π * c / (sqrt(2) * ω) : 2π * c / (nx * ω)
    λz = 2π * c / (sqrt(1 - nx^2) * ω)
    Lx, Lz = nwav * λx, nwav * λz
    lc = λz * nwav / resol
    left_b, interf, right_b = -Lz/2, 0.0, Lz/2
    bot_b, top_b = -Lx/2, Lx/2

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.model.add("box2d")

    # points (z, x) and lines; line 7 is the interface
    gmsh.model.geo.addPoint(left_b,  bot_b, 0, lc, 1)
    gmsh.model.geo.addPoint(interf,  bot_b, 0, lc, 2)
    gmsh.model.geo.addPoint(right_b, bot_b, 0, lc, 3)
    gmsh.model.geo.addPoint(right_b, top_b, 0, lc, 4)
    gmsh.model.geo.addPoint(interf,  top_b, 0, lc, 5)
    gmsh.model.geo.addPoint(left_b,  top_b, 0, lc, 6)
    gmsh.model.geo.addLine(1, 2, 1)
    gmsh.model.geo.addLine(2, 3, 2)
    gmsh.model.geo.addLine(3, 4, 3)
    gmsh.model.geo.addLine(4, 5, 4)
    gmsh.model.geo.addLine(5, 6, 5)
    gmsh.model.geo.addLine(6, 1, 6)
    gmsh.model.geo.addLine(2, 5, 7)

    gmsh.model.geo.addCurveLoop([1, 7, 5, 6], 1)       # medium 1 (z < 0)
    gmsh.model.geo.addPlaneSurface([1], 1)
    gmsh.model.geo.addCurveLoop([2, 3, 4, -7], 2)      # medium 2 (z > 0)
    gmsh.model.geo.addPlaneSurface([2], 2)
    gmsh.model.geo.synchronize()

    gmsh.model.addPhysicalGroup(0, [1, 6], 1)
    gmsh.model.setPhysicalName(0, 1, "LportNodes")
    gmsh.model.addPhysicalGroup(0, [3, 4], 2)
    gmsh.model.setPhysicalName(0, 2, "RportNodes")
    gmsh.model.addPhysicalGroup(1, [6], 3)
    gmsh.model.setPhysicalName(1, 3, "LportEdges")
    gmsh.model.addPhysicalGroup(1, [3], 4)
    gmsh.model.setPhysicalName(1, 4, "RportEdges")
    gmsh.model.addPhysicalGroup(1, [7], 5)
    gmsh.model.setPhysicalName(1, 5, "Interface")
    gmsh.model.addPhysicalGroup(2, [1], 6)
    gmsh.model.setPhysicalName(2, 6, "Medium 1")
    gmsh.model.addPhysicalGroup(2, [2], 7)
    gmsh.model.setPhysicalName(2, 7, "Medium 2")
    gmsh.model.addPhysicalGroup(1, [4, 5], 6)
    gmsh.model.setPhysicalName(1, 6, "top")
    gmsh.model.addPhysicalGroup(1, [1, 2], 7)
    gmsh.model.setPhysicalName(1, 7, "bottom")

    # periodicity: top lines = bottom lines translated by Lx along x
    translation = [1, 0, 0, 0, 0, 1, 0, Lx, 0, 0, 1, 0, 0, 0, 0, 1]
    gmsh.model.mesh.setPeriodic(1, [5], [1], translation)
    gmsh.model.mesh.setPeriodic(1, [4], [2], translation)

    gmsh.model.mesh.generate(2)
    gmsh.write(file)
    gmsh.finalize()
end

if abspath(PROGRAM_FILE) == @__FILE__
    generate_mesh(joinpath(@__DIR__, "box2D_6k.msh"); resol = 50)
    generate_mesh(joinpath(@__DIR__, "box2D_15k.msh"); resol = 81)
end
