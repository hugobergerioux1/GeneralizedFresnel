# Cold-plasma media and their dielectric tensor (time dependence e^{-iωt})

const c₀ = 299792458.0               # speed of light in vacuum [m/s]
const μ₀ = 4π * 1e-7                 # vacuum permeability [H/m]
const ϵ₀ = 8.854187817620389e-12     # vacuum permittivity [F/m]
const Z₀ = sqrt(μ₀ / ϵ₀)             # vacuum impedance [Ω]

"""
    ColdPlasma(; X = 0, Y = 0, ν = 0, Xi = 0, Yi = 0, νi = 0, θB = 0, ϕB = 0)

Homogeneous cold plasma, all frequencies being normalized to the wave frequency ω:

    X  = (ω_pe/ω)²   electron plasma frequency      Xi = (ω_pi/ω)²   ion plasma frequency
    Y  = ω_ce/ω      electron cyclotron frequency   Yi = ω_ci/ω      ion cyclotron frequency
    ν  = ν_e/ω       electron collision frequency   νi = ν_i/ω       ion collision frequency

The static magnetic field B₀ is along (sin θB cos ϕB, sin ϕB, cos θB cos ϕB), the angles being
in degrees: B₀ ∥ z (normal to the interface) for θB = ϕB = 0. `ColdPlasma()` is vacuum.
"""
Base.@kwdef struct ColdPlasma
    X::Float64  = 0.0
    Y::Float64  = 0.0
    ν::Float64  = 0.0
    Xi::Float64 = 0.0
    Yi::Float64 = 0.0
    νi::Float64 = 0.0
    θB::Float64 = 0.0
    ϕB::Float64 = 0.0
end

function Base.show(io::IO, m::ColdPlasma)
    if m.X == 0 && m.Xi == 0
        print(io, "vacuum")
        return
    end
    print(io, "cold plasma, X = $(m.X), Y = $(m.Y)")
    m.ν != 0 && print(io, ", ν = $(m.ν)")
    m.Xi != 0 && print(io, ", Xi = $(m.Xi), Yi = $(m.Yi)", m.νi != 0 ? ", νi = $(m.νi)" : "")
    if m.θB == 0 && m.ϕB == 0
        print(io, " (B₀ ∥ z)")
    else
        print(io, " (θB = $(m.θB)°, ϕB = $(m.ϕB)°)")
    end
    S, _, P = stix_parameters(m)
    if real(S) * real(P) < 0
        @printf(io, ", resonance cones at %.1f° from B₀", atand(sqrt(-real(P) / real(S))))
    end
end

"""
    S, D, P = stix_parameters(medium; Z = 0)

Stix parameters of the medium. `Z` is an additional collisionality ν/ω: it is a
purely numerical damping, used by the FEM solver only.
"""
function stix_parameters(m::ColdPlasma; Z = 0.0)
    ν = m.ν + Z
    Y = abs(m.Y)
    R = 1 - m.X / (1 + im*ν - Y) - m.Xi / (1 + im*m.νi + m.Yi)
    L = 1 - m.X / (1 + im*ν + Y) - m.Xi / (1 + im*m.νi - m.Yi)
    P = 1 - m.X / (1 + im*ν)     - m.Xi / (1 + im*m.νi)
    return (R + L) / 2, (R - L) / 2, P
end

# Rotation R_y(θB) R_x(ϕB) from the frame of B₀ (B₀ ∥ z) to the (x, y, z) frame of the problem
function rotation(m::ColdPlasma)
    θ, ϕ = deg2rad(m.θB), deg2rad(m.ϕB)
    Rθ = [ cos(θ)  0  sin(θ)
           0       1  0
          -sin(θ)  0  cos(θ)]
    Rϕ = [1   0       0
          0   cos(ϕ)  sin(ϕ)
          0  -sin(ϕ)  cos(ϕ)]
    return Rθ * Rϕ
end

# Unit vector along B₀
field_direction(m::ColdPlasma) = rotation(m)[:, 3]

"""
    κ = dielectric_tensor(medium; Z = 0)

Relative dielectric tensor (3×3 complex matrix) in the (x, y, z) frame of the problem: the
Stix tensor, written with B₀ along z, is rotated by R_y(θB) R_x(ϕB).
"""
function dielectric_tensor(m::ColdPlasma; Z = 0.0)
    S, D, P = stix_parameters(m; Z)
    κ_B = [S     -im*D  0
           im*D   S     0
           0      0     P]
    R = rotation(m)
    return R * κ_B * R'
end
