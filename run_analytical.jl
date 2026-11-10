# Analytical solution of the generalized Fresnel problem (no finite elements)
#
#   julia --project=. run_analytical.jl        (or include("run_analytical.jl") from the REPL)

using Pkg
Pkg.activate(@__DIR__)

isdefined(Main, :FresnelAnalytical) || include(joinpath(@__DIR__, "src", "analytical", "FresnelAnalytical.jl"))
using .FresnelAnalytical

case = "R8"                                   # any file of cases/

include(joinpath(@__DIR__, "cases", "$case.jl"))

# the 4 plane waves allowed in each medium (dispersion relation), sorted by energy flow
println("\nWaves in z < 0 ($medium_L):")
print_waves(plane_waves(medium_L, nx, ny))
println("\nWaves in z > 0 ($medium_R):")
print_waves(plane_waves(medium_R, nx, ny))

# amplitudes of the reflected and transmitted waves
sol = solve_fresnel(medium_L, medium_R, nx, ny, A_in)
print_solution(sol)
