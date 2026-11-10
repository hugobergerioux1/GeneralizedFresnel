"""
    FresnelAnalytical

Analytical solution of the generalized Fresnel problem: reflection and transmission of
plane waves at the planar interface z = 0 between two homogeneous media described by the
cold-plasma model (vacuum, isotropic plasma, magnetized plasma, with or without collisions).

    cold_plasma.jl   media and their dielectric tensor
    dispersion.jl    dispersion relation: normal indices nz and polarizations of the waves
    fresnel.jl       attribution of the waves (Poynting flux), interface conditions, report
"""
module FresnelAnalytical

using LinearAlgebra
using Printf

export c₀, μ₀, ϵ₀, Z₀
export ColdPlasma, stix_parameters, dielectric_tensor
export PlaneWave, wave_operator, dispersion_roots, polarization, poynting_z, plane_waves
export solve_fresnel, print_waves, print_solution

include("cold_plasma.jl")
include("dispersion.jl")
include("fresnel.jl")

end
