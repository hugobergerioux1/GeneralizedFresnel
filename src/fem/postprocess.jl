# Post-processing: comparison with the analytical solution, VTU output.

println()
println("Post-processing...")

## Numerical field
if formulation == "SP"
    E_zx, E_y, s_h = xh         # (Ez, Ex), Ey and the Lagrange multiplier (zero for the exact field)
else
    E_zx, E_y = xh
end
Enum = F(E_zx, E_y)             # E = (Ex, Ey, Ez)

## Analytical field: wave i, Asol[i] Pᵢ exp(i kᵢ·r)
in_medium(x, i) = (i <= 4) == (x[1] < 0.0)
wave_phase(x, i) = exp(im*(k[i,1]*x[2] + k[i,3]*x[1]))
Pvec  = [VectorValue(Po[i,1], Po[i,2], Po[i,3]) for i in 1:8]
curlP = [im*VectorValue(cross(k[i,:], Po[i,:])...) for i in 1:8]     # ∇×(P e^{ik·r}) = i k×P e^{ik·r}
E_ana(x)     = sum(in_medium(x,i) ? Asol[i]*Pvec[i]*wave_phase(x,i)  : zero(Pvec[i])  for i in 1:8)
curlE_ana(x) = sum(in_medium(x,i) ? Asol[i]*curlP[i]*wave_phase(x,i) : zero(curlP[i]) for i in 1:8)

## Relative error in the scaled H(curl) norm, the analytical field being evaluated at the quadrature points
Ea = CellField(E_ana, Ω)
curlEa = CellField(curlE_ana, Ω)
sq_norm(x) = norm(x)^2
norm_Hcurl(E, curlE) = sqrt(sum( ∫( sq_norm∘E + (sq_norm∘curlE)/abs2(k₀) )dΩ ))

err_Hcurl = norm_Hcurl(Ea - Enum, curlEa - curl2D(E_zx, E_y)) / norm_Hcurl(Ea, curlEa)

println()
println("-"^60)
@printf("Relative H(curl) error        : %.6e\n", err_Hcurl)
@printf("Assembly + solve time         : %.3f s\n", sim_time)
println("-"^60)

## Outputs
case_name = "$(case)_$(formulation)_$(mesh)$(order)_Z=$(Z)_$(round(Int, nelm/1000))k"

# one line per simulation in a summary table
summary_file = joinpath(outdir, "summary.txt")
isfile(summary_file) || open(io -> println(io, join(["case", "formulation", "mesh", "order", "Z", "nelm",
    "ndofs", "k0_heff", "err_Hcurl", "time_s", "compile_s"], '\t')), summary_file, "w")
open(summary_file, "a") do io
    meshname = mesh == "T" ? basename(meshfile) : "R$(partition[1])x$(partition[2])"
    println(io, join([case, formulation, meshname, order, Z, nelm, ndofs, k₀*heff, err_Hcurl,
                      sim_time, compile_time], '\t'))
end
println("Summary appended to $summary_file")

# fields in ParaView format
if write_vtu
    Eana = F(Gridap.interpolate(x -> (E = E_ana(x); VectorValue(E[3], E[1])), V),
             Gridap.interpolate(x -> E_ana(x)[2], Q))
    fields = Pair{String,Any}[
        "Re(E_num)"     => real(Enum),
        "Im(E_num)"     => imag(Enum),
        "Re(E_ana)"     => real(Eana),
        "Im(E_ana)"     => imag(Eana),
        "|E_num|"       => MAG(Enum),           # component-wise modulus
        "|E_ana|"       => MAG(Eana),
        "arg(E_num)"    => phase(Enum),         # component-wise phase
        "arg(E_ana)"    => phase(Eana),
        "|E_num-E_ana|" => MAG(Enum - Eana),
    ]
    formulation == "SP" && push!(fields, "|s|" => Operation(abs)(s_h))
    writevtk(Ω, joinpath(outdir, case_name * ".vtu"), cellfields=fields)
    println("Fields written to $(joinpath(outdir, case_name * ".vtu"))")
end

