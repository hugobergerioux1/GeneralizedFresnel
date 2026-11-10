# Finite element simulation of the generalized Fresnel problem, compared with the analytical
# solution. Included by run_fem.jl, which defines the parameters
#   case, formulation, mesh, meshfile, partition, order, Z, nthreads, outdir, write_vtu

using Gridap
using GridapGmsh
using Gridap.Fields
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.ReferenceFEs
using Gridap.Arrays
using Gridap.CellData
using LinearAlgebra
using SparseArrays
using Printf

BLAS.set_num_threads(nthreads)

isdefined(Main, :FresnelAnalytical) || include(joinpath(@__DIR__, "..", "analytical", "FresnelAnalytical.jl"))
using .FresnelAnalytical

include(joinpath(@__DIR__, "operators.jl"))

## 1. Analytical solution: wave vectors, polarizations and amplitudes of the 8 plane waves.
include(joinpath(@__DIR__, "..", "..", "cases", "$case.jl"))
ny == 0 || error("the 2D FE model assumes fields independent of y: set ny = 0 in cases/$case.jl")
sol = solve_fresnel(medium_L, medium_R, nx, ny, A_in)
print_solution(sol)

println("\nNUMERICAL SOLUTION    case $case, formulation $formulation, mesh $mesh, order $order")

ω  = 2π * freq                  # angular frequency [rad/s]
k₀ = ω / c₀                     # vacuum wavenumber [1/m]
Asol = sol.A                    # amplitudes of the 8 waves (at z = 0)
Po   = sol.P                    # polarizations, Po[i, :] = (Px, Py, Pz) of wave i
k    = k₀ * sol.n               # wave vectors [1/m], k[i, :] = (kx, ky, kz) of wave i

# Dielectric tensor of the FE model: the numerical collisionality Z is added to the media
κ_numL = dielectric_tensor(medium_L; Z)
κ_numR = dielectric_tensor(medium_R; Z)
function κn(x)
    if x[1] <= 0
        return TensorValue(κ_numL)
    else
        return TensorValue(κ_numR)
    end
end

## 2. Mesh and mixed FE spaces
include(joinpath(@__DIR__, "discretization.jl"))

## 3. Waveport boundary conditions
include(joinpath(@__DIR__, "waveports.jl"))

## 4. Assembly and solution (timed)
include(joinpath(@__DIR__, "solve.jl"))

## 5. Errors, VTU output
include(joinpath(@__DIR__, "postprocess.jl"))
