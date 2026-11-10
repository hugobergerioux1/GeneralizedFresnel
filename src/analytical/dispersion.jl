# Dispersion relation of a homogeneous medium: for given tangential indices (nx, ny), the 4
# plane waves E e^{i(k₀ n·r − ωt)}, n = (nx, ny, nz), allowed by Maxwell's equations.

"""
    PlaneWave(n, E)

Plane wave E e^{i k₀ n·r}: refractive index vector n = k/k₀ (nz is complex for evanescent
or damped waves) and polarization vector E.
"""
struct PlaneWave
    n::Vector{ComplexF64}
    E::Vector{ComplexF64}
end

"""
    M = wave_operator(κ, n)

M(n) = κ + n nᵀ − (n·n) I. The plane wave of index n and polarization E satisfies Maxwell's
equations in the medium κ iff M(n) E = 0
"""
wave_operator(κ, n) = κ + n * transpose(n) - (transpose(n) * n) * I

"""
    nz = dispersion_roots(κ, nx, ny)

The 4 roots nz of the dispersion relation det M(nx, ny, nz) = 0 (kz = k₀ nz).
"""
function dispersion_roots(κ, nx, ny)
    # det M is a polynomial of degree 4 in nz: its coefficients are recovered
    # exactly from 5 samples on the unit circle (discrete Fourier transform)
    ζ = [cispi(2j / 5) for j in 0:4]
    d = [det(wave_operator(κ, [nx, ny, z])) for z in ζ]
    c = [sum(d .* conj.(ζ) .^ j) / 5 for j in 0:4]          # det M = Σⱼ c[j+1] nz^j
    abs(c[5]) > 1e-12 * maximum(abs, c) ||
        error("κ_zz = 0: the dispersion relation is degenerate (cutoff P = 0 with B₀ ∥ z)")
    # the roots are the eigenvalues of the companion matrix
    C = diagm(-1 => ones(ComplexF64, 3))
    C[:, 4] = -c[1:4] / c[5]
    return eigvals(C)
end

"""
    E = polarization(κ, n)

Polarization of the wave of index n: unit null vector of M(n) (right singular vector of the
smallest singular value), with its largest component made real and positive.
"""
function polarization(κ, n)
    E = svd(wave_operator(κ, n)).V[:, end]
    j = argmax(abs.(E))
    return E * (abs(E[j]) / E[j])
end

"""
    Sz = poynting_z(waves, A)
    Sz = poynting_z(wave)

z-component at z = 0 of the time-averaged Poynting vector
"""
function poynting_z(waves, A)
    E = sum(A[i] * waves[i].E for i in eachindex(waves))
    H = sum(A[i] * cross(waves[i].n, waves[i].E) for i in eachindex(waves))
    return real(cross(E, conj(H))[3])
end
poynting_z(w::PlaneWave) = poynting_z([w], [1])

# Direction (±1) in which the wave carries energy along z: an evanescent or damped wave
# decays in the direction of its energy flow, a propagating wave follows its Poynting flux.
direction(w::PlaneWave) = imag(w.n[3]) != 0 ? sign(imag(w.n[3])) : sign(poynting_z(w))

# Rotation sense of the polarization about the unit vector b: i(E × E*)·b / |E|², from +1
# (right-handed, i.e. rotating like the electrons about B₀) to −1 (left-handed), 0 if linear.
handedness(E, b) = real(im * dot(b, cross(E, conj(E)))) / norm(E)^2

"""
    waves = plane_waves(medium, nx, ny)
    waves = plane_waves(κ, nx, ny; b = [0, 0, 1])

The 4 plane waves of the medium (or of the dielectric tensor κ, B₀ being along b) for the
tangential indices (nx, ny), ordered as
    [forward 1, forward 2, backward 1, backward 2]

Within a pair, the R-like wave comes first and the L-like one second; waves with the same 
rotation sense are ordered by decreasing |Re nz|. In isotropic media the waves of a pair 
are degenerate: TE first, then TM.
"""
plane_waves(m::ColdPlasma, nx, ny) = plane_waves(dielectric_tensor(m), nx, ny; b = field_direction(m))

function plane_waves(κ::AbstractMatrix, nx, ny; b = [0.0, 0.0, 1.0], tol = 1e-9)
    isotropic = norm(κ - κ[1, 1] * I) <= tol * norm(κ)
    if isotropic
        waves = isotropic_waves(κ[1, 1], nx, ny)
    else
        waves = PlaneWave[]
        for nz in dispersion_roots(κ, nx, ny)
            abs(imag(nz)) < tol * max(1, abs(nz)) && (nz = complex(real(nz)))   # round-off
            n = [nx, ny, nz]
            push!(waves, PlaneWave(n, polarization(κ, n)))
        end
    end
    fwd = filter(w -> direction(w) > 0, waves)
    bwd = filter(w -> direction(w) < 0, waves)
    length(fwd) == length(bwd) == 2 ||
        error("the 4 waves cannot be split into 2 forward and 2 backward waves")
    if !isotropic
        precedes(fwd[1], fwd[2], b) || reverse!(fwd)
        precedes(bwd[1], bwd[2], b) || reverse!(bwd)
    end
    return [fwd; bwd]
end

# Order within a pair: R-like before L-like, otherwise larger |Re nz| first
function precedes(w1, w2, b)
    h1, h2 = handedness(w1.E, b), handedness(w2.E, b)
    abs(h1 - h2) > 1e-6 && return h1 > h2
    return abs(real(w1.n[3])) > abs(real(w2.n[3]))
end

"""
    waves = isotropic_waves(ε, nx, ny)

Waves of an isotropic medium κ = ε I. Both roots nz = ±√(ε − nx² − ny²) are double and carry
a TE wave, E = n̂ₜ × ẑ (i.e. (0, −1, 0) for ny = 0), and a TM wave, E ∝ σ (nz n̂ₜ − nₜ ẑ),
where n̂ₜ is the unit tangential direction and σ = ±1 the direction of propagation: the TM
waves of a pair are thus mirror images of each other through the interface.
"""
function isotropic_waves(ε, nx, ny)
    nt = hypot(nx, ny)
    t = nt > 0 ? [nx, ny, 0] / nt : [1.0, 0.0, 0.0]
    E_TE = complex([t[2], -t[1], 0.0])
    s = sqrt(complex(ε - nt^2))
    waves = PlaneWave[]
    for nz in (s, -s)
        n = [nx, ny, nz]
        σ = direction(PlaneWave(n, E_TE))
        E_TM = σ * (nz * t - [0, 0, nt])
        push!(waves, PlaneWave(n, E_TE), PlaneWave(n, E_TM / norm(E_TM)))
    end
    return waves
end

"""
    print_waves(waves)

Prints plane waves: index, polarization, Poynting flux.
"""
function print_waves(waves)
    for (i, w) in enumerate(waves)
        @printf("  %d   nz = %-18s P = %-48s Sz = %+.4e   energy towards %sz\n", i,
                cstr(w.n[3]), vstr(w.E), poynting_z(w), direction(w) > 0 ? "+" : "-")
    end
end

# compact strings of complex numbers and vectors, negligible parts dropped
function cstr(z; d = 4)
    tol = 0.5 * 10.0^(-d)
    re = abs(real(z)) < tol ? 0.0 : real(z)
    im_ = abs(imag(z)) < tol ? 0.0 : imag(z)
    im_ == 0 && return @sprintf("%.*f", d, re)
    re == 0 && return @sprintf("%.*fi", d, im_)
    return @sprintf("%.*f%s%.*fi", d, re, im_ < 0 ? "-" : "+", d, abs(im_))
end
vstr(v; d = 4) = "(" * join(cstr.(v; d), ", ") * ")"
