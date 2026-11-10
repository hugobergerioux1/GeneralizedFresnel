# Generalized Fresnel problem: interface z = 0 between medium_L (z < 0) and medium_R (z > 0).

const WAVE_ROLES = ("incoming", "incoming", "reflected", "reflected",
                    "incoming", "incoming", "transmitted", "transmitted")

"""
    c = interface_conditions(wave, κ)

Quantities that must be continuous across the interface z = 0:

    c = (Ex, Ey, Dz/ϵ₀, cBz, cBx, cBy),   with D = ϵ₀ κ E and c B = n × E,

i.e. the tangential components of E, the normal component of D and the magnetic field B.
"""
function interface_conditions(w, κ)
    cB = cross(w.n, w.E)
    return [w.E[1], w.E[2], (κ * w.E)[3], cB[3], cB[1], cB[2]]
end

"""
    sol = solve_fresnel(medium_L, medium_R, nx, ny, A_in)

Plane waves of tangential indices (nx, ny) at the interface z = 0 between medium_L (z < 0)
and medium_R (z > 0). The 8 waves are numbered

            z < 0                              z > 0
    1, 2  incoming  (energy → +z)      5, 6  incoming    (energy → −z)
    3, 4  reflected (energy → −z)      7, 8  transmitted (energy → +z)

(see `plane_waves` for the attribution and the order within each pair), the field of wave i
being Aᵢ Pᵢ e^{i k₀ nᵢ·r}. Given the incoming amplitudes A_in = [A1, A2, A5, A6], the outgoing
ones (A3, A4, A7, A8) follow from the interface conditions at z = 0.

Returns a named tuple with the 8 `waves`, their amplitudes `A`, the index vectors `n` and
polarizations `P` as 8×3 matrices (row i = wave i), the dielectric tensors `κ_L` and `κ_R`,
the incoming power flux `S_in` and the reflected and transmitted fractions `R` and `T` of the incoming power.
"""
function solve_fresnel(medium_L, medium_R, nx, ny, A_in)
    κ_L = dielectric_tensor(medium_L)
    κ_R = dielectric_tensor(medium_R)
    wL = plane_waves(medium_L, nx, ny)                
    wR = plane_waves(medium_R, nx, ny)
    waves = [wL[1], wL[2], wL[3], wL[4], wR[3], wR[4], wR[1], wR[2]]

    # interface conditions of each wave (waves 1-4 in medium_L, waves 5-8 in medium_R)
    c = [interface_conditions(waves[i], i <= 4 ? κ_L : κ_R) for i in 1:8]

    # continuity at z = 0 with the unknown outgoing amplitudes on the left 
    # and the given incoming ones on the right:
    #   M_out [A3, A4, A7, A8] = M_in [A1, A2, A5, A6]
    M_out = hcat(-c[3], -c[4], c[7], c[8])
    M_in  = hcat( c[1],  c[2], -c[5], -c[6])

    # 6 conditions for 4 unknowns: for waves satisfying Maxwell's equations, the continuity of
    # Dz and Bz follows from that of the tangential fields
    A_out = M_out \ (M_in * A_in)
    # residual printed as a check
    residual = norm(M_out * A_out - M_in * A_in) / norm(M_in * A_in)

    A = zeros(ComplexF64, 8)
    A[[1, 2, 5, 6]] = A_in
    A[[3, 4, 7, 8]] = A_out

    # power balance at the interface (R + T = 1 for lossless media)
    S_in = poynting_z(waves[1:2], A[1:2]) - poynting_z(waves[5:6], A[5:6])
    R = -poynting_z(waves[3:4], A[3:4]) / S_in
    T = poynting_z(waves[7:8], A[7:8]) / S_in

    n = permutedims(reduce(hcat, [w.n for w in waves]))
    P = permutedims(reduce(hcat, [w.E for w in waves]))
    return (; medium_L, medium_R, nx, ny, κ_L, κ_R, waves, A, n, P, S_in, R, T, residual)
end

"""
    print_solution(sol)

Prints the 8 waves of a solution of `solve_fresnel` (index, polarization, amplitude, power
flux) and the power balance at the interface.
"""
function print_solution(sol)
    println("\n" * "="^100)
    println("ANALYTICAL SOLUTION     field of wave i: Aᵢ Pᵢ exp(i k₀ nᵢ·r),  nᵢ = (nx, ny, nzᵢ)")
    println("="^100)
    println("  z < 0 : ", sol.medium_L)
    println("  z > 0 : ", sol.medium_R)
    @printf("  nx = %.6g, ny = %.6g\n\n", sol.nx, sol.ny)
    @printf("  %-4s %-6s %-12s %-18s %-48s %-18s %s\n",
            "wave", "side", "role", "nz", "P = (Px, Py, Pz)", "A", "|A|²Sz/S_in")
    for i in 1:8
        w = sol.waves[i]
        @printf("  %-4d %-6s %-12s %-18s %-48s %-18s %+.6f\n", i, i <= 4 ? "z < 0" : "z > 0",
                WAVE_ROLES[i], cstr(w.n[3]), vstr(w.E), cstr(sol.A[i]),
                abs2(sol.A[i]) * poynting_z(w) / sol.S_in)
    end
    @printf("\n  power fractions:  reflected R = %.6f,  transmitted T = %.6f,  R + T = %.6f",
            sol.R, sol.T, sol.R + sol.T)
    println("  (= 1 for lossless media)")
    @printf("  interface conditions (Ex, Ey, Dz, Bz, Bx, By): relative residual %.1e\n",
            sol.residual)
    println("="^100)
end
