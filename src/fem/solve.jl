# Weak form of ∇×∇×E − k₀² κ E = 0 (fields independent of y) with the waveport terms,
# assembled and solved by sparse LU. This block is timed: the time reported is that of the
# assembly and solution, without Julia's just-in-time compilation.
#
#   W:  standard formulation of the wave equation
#   SP: saddle-point formulation, the constraint ∇·(κE) = 0 being enforced weakly with the
#       Lagrange multiplier s

println()
println("Building general weak form...")
timing = @timed begin
    if formulation == "SP"
        println("   with saddle-point formulation...")
        a((u,p,s), (v,q,t)) = ∫( curl2D(v,q)'⋅curl2D(u,p) - (k₀^2)*F(v,q)'⋅((κn∘(Px))⋅F(u,p))
            -  F(v,q)'⋅((κn∘(Px))⋅grad2D(s)) + grad2D(t)'⋅((κn∘(Px))⋅F(u,p)) )dΩ
        b((v,q,t)) = ∫(im*0*q)dΩ
        op = AffineFEOperator(a,b,X,Y)
    else
        println("   with wave equation formulation...")
        a((u,p), (v,q)) = ∫( curl2D(v,q)'⋅curl2D(u,p) - (k₀^2)*F(v,q)'⋅((κn∘(Px))⋅F(u,p)) )dΩ
        b((v,q)) = ∫(im*0*q)dΩ
        op = AffineFEOperator(a,b,X,Y)
    end
    op.op.matrix .+= PEᵣᵢ - PEₗₑ      # waveports
    op.op.vector .+= Uᵣᵢ - Uₗₑ

    println()
    println("Solving...")
    xh = Gridap.solve(op)
end

sim_time = timing.time - timing.compile_time    
compile_time = timing.compile_time
@printf("Assembly + solve time: %.3f s\n", sim_time)
