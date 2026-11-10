# GeneralizedFresnel

**Reflection and transmission of electromagnetic plane waves at a planar interface between two
cold-plasma media: analytical solver and 2D finite element solver, in Julia.**

This repository contains the code used for the simulations of the article

> H. Bergerioux and M. Merino, *On Spurious Mesh-Scale Oscillations Associated with
> Resonance-Cone Regimes in Wave-Plasma Simulations*, arXiv:2610.31021 (2026).
> [https://arxiv.org/abs/2610.31021](https://arxiv.org/abs/2610.31021)

where the problem, the methods and the results are described in detail. The code is made of two
blocks, which can be used independently:

- an **analytical solver** (`src/analytical/`, plain Julia), for two homogeneous cold-plasma
  media (vacuum, isotropic or magnetized plasma, with or without collisions): waves allowed in
  each medium, attribution to the incoming or outgoing directions from their energy flux, and
  amplitudes of the reflected and transmitted waves;
- a **2D finite element solver** (`src/fem/`, based on [Gridap.jl](https://github.com/gridap/Gridap.jl))
  of the same problem: mixed Nédélec–Lagrange elements on triangles or rectangles, standard (W)
  or saddle-point (SP) formulation, periodic lateral boundaries and waveport boundary conditions
  built from the analytical waves. The numerical solution is compared with the analytical one.

## Repository structure

```
GeneralizedFresnel/
├── run_analytical.jl          analytical solution only
├── run_fem.jl                 finite element simulation, compared with the analytical solution
├── cases/                     test cases: media, frequency, incidence, incoming amplitudes
├── meshes/
│   ├── box2D_6k.msh           triangular meshes of the article (Gmsh)
│   ├── box2D_15k.msh
│   └── generate_mesh.jl       mesh generator
├── src/
│   ├── analytical/            module FresnelAnalytical
│   │   ├── FresnelAnalytical.jl  module definition
│   │   ├── cold_plasma.jl     media, cold-plasma dielectric tensor
│   │   ├── dispersion.jl      dispersion relation: waves (nz, polarization) of a medium
│   │   └── fresnel.jl         wave attribution, interface conditions, printed report
│   └── fem/
│       ├── simulation.jl      main script of a simulation, runs the files below in order
│       ├── operators.jl       2D curl and gradient, 2D → 3D field reconstruction
│       ├── discretization.jl  mesh, integration measures, mixed FE spaces
│       ├── waveports.jl       waveport boundary conditions
│       ├── solve.jl           weak form, assembly and solution (timed)
│       └── postprocess.jl     error, VTU output
├── Project.toml               dependencies
└── Manifest.toml              exact versions of the dependencies used for the article
```

## Installation

The code requires [Julia](https://julialang.org/downloads/) 1.12 (the results of the article were
obtained with Julia 1.12.6). Clone the repository and install the dependencies:

```bash
git clone https://github.com/hugobergerioux1/GeneralizedFresnel.git
cd GeneralizedFresnel
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

`Manifest.toml` pins the exact versions of the packages used for the article (Gridap 0.19.8,
GridapGmsh 0.7.3, which also provides Gmsh 4.15).

## Quick start

```bash
julia --project=. run_analytical.jl     # analytical solution
julia --project=. run_fem.jl            # finite element simulation (results in results/)
```

The case and the other parameters are set at the top of each script. The scripts can also be
run from a Julia REPL started with `julia --project=.`, with `include("run_fem.jl")`: the
following runs in the same session are much faster, since the code is already compiled.

## Conventions

The interface $z = 0$ separates `medium_L` ($z < 0$) from `medium_R` ($z > 0$), and the fields
vary in time as $e^{-i\omega t}$. All the waves share the tangential refractive index
$(n_x, n_y)$, and wave $i$ is the plane wave $A_i\, \mathbf P_i\, e^{i k_0 \mathbf n_i\cdot\mathbf r}$,
with $\mathbf n_i = (n_x, n_y, n_{z,i})$, amplitude $A_i$ (at $z = 0$) and unit polarization
$\mathbf P_i$. The 8 waves are numbered as follows:

|          | $z < 0$                           | $z > 0$                             |
|----------|-----------------------------------|-------------------------------------|
| incoming | 1, 2 (energy towards $+z$)        | 5, 6 (energy towards $-z$)          |
| outgoing | 3, 4 reflected (towards $-z$)     | 7, 8 transmitted (towards $+z$)     |

The incoming amplitudes `A_in = [A1, A2, A5, A6]` are data of the problem. With vacuum on the
left, A1 and A2 are the TE ($\mathbf P = (0, -1, 0)$) and TM components of the incident wave.

A medium is defined by `ColdPlasma(; X, Y, ν, Xi, Yi, νi, θB, ϕB)`, all frequencies being
normalized to the wave frequency $\omega$ (all parameters default to 0; `ColdPlasma()` is vacuum):

| parameter    | meaning                                                              |
|--------------|----------------------------------------------------------------------|
| `X`, `Xi`    | $(\omega_{pe}/\omega)^2$, $(\omega_{pi}/\omega)^2$: electron and ion plasma frequencies |
| `Y`, `Yi`    | $\omega_{ce}/\omega$, $\omega_{ci}/\omega$: electron and ion cyclotron frequencies    |
| `ν`, `νi`    | $\nu_e/\omega$, $\nu_i/\omega$: electron and ion collision frequencies               |
| `θB`, `ϕB`   | orientation of $\mathbf B_0$ [deg], $\mathbf B_0 \parallel (\sin\theta_B\cos\phi_B, \sin\phi_B, \cos\theta_B\cos\phi_B)$ |

## Analytical solver

`run_analytical.jl` prints the 4 waves of each medium ($n_z$, polarization, Poynting flux) and
the solution: amplitudes of the 8 waves, reflected and transmitted power fractions R and T
(R + T = 1 for lossless media), and the residual of the interface conditions.

- **Dispersion relation.** For given $(n_x, n_y)$, the 4 roots $n_z$ are obtained as the
  eigenvalues of the companion matrix of the quartic dispersion polynomial: no initial guess is
  needed, and propagating, evanescent and damped waves are all found. The polarizations are
  computed by SVD (`dispersion.jl`).
- **Attribution of the waves.** A wave goes towards $+z$ if it decays towards $+z$, or, when it
  propagates, if its Poynting flux $S_z$ is positive. Backward waves are thus correctly
  attributed (e.g. the transmitted wave $n_z = -3.12$ of case R3).
- **Order within a pair.** In a magnetized plasma, the R-like wave (polarization rotating like the
  electrons about $\mathbf B_0$) comes first; in vacuum, the TE wave comes first, then the TM wave.
  The polarizations of plasma waves are normalized with their largest component real and positive.
- **Interface conditions.** The continuity of the tangential $\mathbf E$, of the normal
  $\mathbf D$ and of $\mathbf B$ gives the system $M_\text{out} A_\text{out} = M_\text{in} A_\text{in}$
  (`fresnel.jl`).

The solver can also be used from Julia:

```julia
include("src/analytical/FresnelAnalytical.jl")
using .FresnelAnalytical

plasma = ColdPlasma(X = 8.0, Y = 4.0)                      # B₀ ∥ z
waves = plane_waves(plasma, sind(45), 0.0)                # waves[i].n, waves[i].E
sol = solve_fresnel(ColdPlasma(), plasma, sind(45), 0.0, [1.0, 1.0, 0.0, 0.0])
sol.A, sol.R, sol.T                                        # amplitudes, power fractions
print_solution(sol)
```

## Finite element solver

The fields are independent of $y$ ($n_y = 0$) and are computed in the $(z, x)$ plane:
$(E_z, E_x)$ with Nédélec elements of order `order − 1`, $E_y$ with Lagrange elements of order
`order`. `run_fem.jl` first runs the analytical solver, which provides the reference solution and
the modes of the waveports, then:

- **Formulation** (`solve.jl`): `"W"`, standard formulation of the wave equation, or `"SP"`,
  saddle-point formulation where the constraint $\nabla\cdot(\kappa\mathbf E) = 0$ is enforced
  with a Lagrange multiplier.
- **Boundary conditions**: periodic on the top and bottom sides; waveports on the left and right
  sides (`waveports.jl`), which impose the incoming waves and are transparent to the outgoing
  ones, so that the bounded domain replicates the two semi-infinite media.
- **Meshes**: `"T"`, triangles read from a Gmsh file; `"R"`, rectangles of a structured grid
  built by Gridap (see [Meshes](#meshes)).
- **Numerical collisionality** `Z`: artificial electron collision frequency $\nu/\omega$, added
  to the media of the finite element model only (`Z = 0`: exact, lossless problem).

### Parameters

| parameter     | default                  | description                                                         |
|---------------|--------------------------|---------------------------------------------------------------------|
| `case`        | `"R8"`                   | test case, file of `cases/`                                         |
| `formulation` | `"W"`                    | `"W"`: wave equation, `"SP"`: saddle point                          |
| `mesh`        | `"T"`                    | `"T"`: triangles (Gmsh file), `"R"`: rectangles (structured grid)   |
| `meshfile`    | `meshes/box2D_15k.msh`   | Gmsh mesh for `mesh = "T"`                                          |
| `partition`   | `(124, 124)`             | cells along $(z, x)$ for `mesh = "R"`                               |
| `order`       | `2`                      | order of the Lagrange elements (Nédélec: `order − 1`)               |
| `Z`           | `0.0`                    | numerical collisionality $\nu/\omega$                               |
| `nthreads`    | `10`                     | number of BLAS threads (10 for the article)                         |
| `outdir`      | `results/`               | output folder                                                       |
| `write_vtu`   | `true`                   | write the fields in ParaView format                                 |

### Outputs

- **Terminal**: analytical solution, mesh, number of degrees of freedom, relative error and
  computational time.
- `results/summary.txt`: one line per simulation (appended), with the columns `case`,
  `formulation`, `mesh`, `order`, `Z`, `nelm`, `ndofs`, `k0_heff`, `err_Hcurl`, `time_s` and
  `compile_s`.
- `results/<case>_<formulation>_<mesh><order>_Z=<Z>_<N>k.vtu` (ParaView): real and imaginary
  parts, modulus and phase of the numerical and analytical fields, modulus of their difference
  and, for SP, modulus of the multiplier.
- `results/mesh_<mesh>_<N>k_*.vtu`: the mesh and its tags.

**Error.** `err_Hcurl` is the relative error in the $H(\mathrm{curl})$ norm scaled by $k_0$,
$\lVert \mathbf E \rVert_{H(\mathrm{curl})}^2 = \lVert \mathbf E \rVert^2 + k_0^{-2}\lVert \nabla\times\mathbf E \rVert^2$,
the analytical field being evaluated at the quadrature points.

**Computational time.** `time_s` is the time of the timed block of `solve.jl` (assembly,
waveport terms and sparse LU solution), without Julia's just-in-time compilation, which is given
separately in `compile_s`.

## Test cases

All the cases are at 5.8 GHz, with a TE + TM wave (`A_in = [1, 1, 0, 0]`) incident from vacuum
at 45°, and a lossless electron plasma with $\mathbf B_0$ normal to the interface. The cases are
named after the region of the CMA diagram of the plasma; R and T are given by the analytical
solver.

| case     | $z > 0$          | X   | Y   | R        | T        |
|----------|------------------|-----|-----|----------|----------|
| `vacuum` | vacuum           | 0   | 0   | 0        | 1        |
| `R1`     | CMA region I     | 0.3 | 0.3 | 0.043162 | 0.956838 |
| `R3`     | CMA region III   | 0.9 | 0.9 | 0.535179 | 0.464821 |
| `R7`     | CMA region VII   | 2.0 | 4.0 | 0.068952 | 0.931048 |
| `R8`     | CMA region VIII  | 8.0 | 4.0 | 0.370379 | 0.629621 |

A new case is a new file in `cases/` defining `freq` (Hz), `medium_L`, `medium_R`, `nx`, `ny`
and `A_in`. The finite element solver requires `ny = 0` and a mesh consistent with `freq` and
`nx`.

**Reference results** of `run_fem.jl`, case R8, `order = 2`, `Z = 0`:

| formulation | mesh                        | cells | DOFs    | `err_Hcurl` |
|-------------|-----------------------------|-------|---------|-------------|
| W           | T, `box2D_15k.msh`          | 15356 | 107 816 | 4.2049e-2   |
| SP          | T, `box2D_15k.msh`          | 15356 | 138 366 | 4.2049e-2   |
| W           | R, `partition = (124, 124)` | 15376 | 185 008 | 8.3279e-3   |
| SP          | R, `partition = (124, 124)` | 15376 | 246 264 | 8.3280e-3   |

## Meshes

The domain spans 5 wavelengths of the incident wave in each direction, with the interface
$z = 0$ in the middle, and is periodic in $x$. A mesh is therefore valid for one frequency and
one incidence only.

- **Triangles** (`mesh = "T"`): `meshes/box2D_6k.msh` (5850 triangles) and
  `meshes/box2D_15k.msh` (15356 triangles) are the meshes of the article (5.8 GHz, 45°). Other
  meshes are generated with `meshes/generate_mesh.jl`, where `resol` sets the element size:

  ```julia
  include("meshes/generate_mesh.jl")
  generate_mesh("meshes/my_mesh.msh"; freq = 5.8e9, nx = sind(30), nwav = 5, resol = 81)
  ```

  Any Gmsh mesh can be used, provided that the interface is the line $z = 0$, that the top and
  bottom sides are periodic, and that the ports are tagged `LportEdges`, `LportNodes` (left)
  and `RportEdges`, `RportNodes` (right).
- **Rectangles** (`mesh = "R"`): no file needed, the grid is built from `partition`, for the
  frequency and incidence of the case. The number of cells along $z$ must be even for the
  interface to be a grid line.

## Citation

If you use this code, please cite:

```bibtex
@misc{bergerioux2026spurious,
  title         = {On Spurious Mesh-Scale Oscillations Associated with Resonance-Cone Regimes in Wave-Plasma Simulations},
  author        = {Bergerioux, Hugo and Merino, Mario},
  year          = {2026},
  eprint        = {2610.31021},
  archivePrefix = {arXiv},
  primaryClass  = {physics.plasm-ph},
  doi           = {10.48550/arXiv.2610.31021},
  url           = {https://arxiv.org/abs/2610.31021}
}
```

## License

This code is released under the MIT License, see [LICENSE](LICENSE).
