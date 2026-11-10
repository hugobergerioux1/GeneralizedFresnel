# Finite element solution of the generalized Fresnel problem, compared with the analytical one
#
#   julia --project=. run_fem.jl        (or include("run_fem.jl") from the REPL)

using Pkg
Pkg.activate(@__DIR__)

case        = "R8"                                            # any file of cases/
formulation = "W"                                             # "W": wave equation, "SP": saddle point
mesh        = "T"                                             # "T": triangles, "R": rectangles
meshfile    = joinpath(@__DIR__, "meshes", "box2D_15k.msh")   # "T": Gmsh mesh, must match freq and nx
partition   = (124, 124)                                      # "R": cells along (z, x) of a 5λz × 5λx grid
order       = 2                                               # Lagrange order (Nédélec: order - 1)
Z           = 0.0                                             # numerical collisionality ν/ω (FEM only)
nthreads    = 10                                              # BLAS threads (10 for the times of the paper)
outdir      = joinpath(@__DIR__, "results")
write_vtu   = true                                            # fields for ParaView

include(joinpath(@__DIR__, "src", "fem", "simulation.jl"))
