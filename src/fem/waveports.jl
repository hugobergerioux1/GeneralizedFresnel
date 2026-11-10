# Waveport boundary conditions on the left (Γₗ, z = left_b) and right (Γᵣ, z = right_b) ports.
#
# On a port, the field is expanded on the plane waves (modes) of the adjacent medium given by
# the analytical solver: the 2 incoming waves, of known amplitudes, and the 2 outgoing ones,
# whose amplitudes are obtained by projecting the FE field on their tangential traces. The
# boundary term of the weak form, ∫ (n × ∇×E)·v̄ dΓ, is then written with these modes: the
# outgoing waves give the matrix terms PE, the incoming ones the right-hand side terms U.
# The ports are thus transparent to the outgoing waves, which leave the domain without
# reflection: the bounded domain replicates the two semi-infinite media of the analytical
# problem.

# Tangential trace of wave i on a port, for the unknowns ((Ez, Ex), Ey) and without the factor
# exp(i kz z) (constant on a port)...
function modeN(x,i)
    return VectorValue(Po[i,3],Po[i,1])*exp(im*k[i,1]*x[2])
end
function modeL(x,i)
    return Po[i,2]*exp(im*k[i,1]*x[2])
end
# ... and trace of (−ẑ) × ∇×E
function fN(x,i)
    return VectorValue(0.0*im,im*(k[i,3]*Po[i,1]-k[i,1]*Po[i,3]))*exp(im*k[i,1]*x[2])
end
function fL(x,i)
    return im*k[i,3]*Po[i,2]*exp(im*k[i,1]*x[2])
end

println()
println("Computing modes dofs...")

modesNe = Dict{Int, Function}()
for i = 1:8
    modesNe[i] = x -> modeN(x,i)
end
modesLa = Dict{Int, Function}()
for i = 1:8
    modesLa[i] = x -> modeL(x,i)
end

modofs = SparseVector{ComplexF64, Int64}[]
if formulation == "SP"          # the multiplier component of the modes is 0
    for i = 1:8
        push!(modofs, sparse(Gridap.interpolate([modesNe[i],modesLa[i],0],Y).free_values))
    end
else
    for i = 1:8
        push!(modofs, sparse(Gridap.interpolate([modesNe[i],modesLa[i]],Y).free_values))
    end
end

fNe = Dict{Int, Function}()
for i = 1:8
    fNe[i] = x -> fN(x,i)
end
fLa = Dict{Int, Function}()
for i = 1:8
    fLa[i] = x -> fL(x,i)
end

fdofs = SparseVector{ComplexF64, Int64}[]
if formulation == "SP"
    for i=1:8
        push!(fdofs, sparse(Gridap.interpolate([fNe[i],fLa[i],0],Y).free_values))
    end
else
    for i=1:8
        push!(fdofs, sparse(Gridap.interpolate([fNe[i],fLa[i]],Y).free_values))
    end
end

# projection on the tangential components (Ex, Ey)
TΓ = TensorValue(
    1.0+0im, 0.0+0im, 0.0+0im,
    0.0+0im, 1.0+0im, 0.0+0im,
    0.0+0im, 0.0+0im, 0.0+0im
)

# amplitudes of the incoming waves at the ports
A_inL = [Asol[i]*exp(im*k[i,3]*left_b) for i in 1:2]      # waves 1, 2 at z = left_b
A_inR = [Asol[i]*exp(im*k[i,3]*right_b) for i in 5:6]     # waves 5, 6 at z = right_b

## LEFT PORT: incoming waves 1, 2 and outgoing (reflected) waves 3, 4
println("Building left port...")

if formulation == "SP"
    aₗₑ((u,p,s),(v,q,t)) = ∫( F(v,q)'⋅ (TΓ ⋅ F(u,p)) )dΓₗ
else
    aₗₑ((u,p),(v,q)) = ∫( F(v,q)'⋅ (TΓ ⋅ F(u,p)) )dΓₗ
end
Mₗₑ = assemble_matrix(aₗₑ,X,Y)

B_L = zeros(Complex{Float64}, 2, 2)      # projections of the outgoing modes on themselves
A_L = zeros(Complex{Float64}, 2, 2)      # and on the incoming modes
αL = zeros(Complex{Float64}, 2, 2)
for i = 1:2
    for m = 1:2
        B_L[i,m] = modofs[i+2]'*Mₗₑ*modofs[m+2]
        A_L[i,m] = modofs[i+2]'*Mₗₑ*modofs[m]
    end
end
BinvL = inv(B_L)

for i = 1:2
    for m = 1:2
        αL[i,m] = BinvL[i,m]*(A_L[m,1]*A_inL[1] + A_L[m,2]*A_inL[2])
    end
end

PEₗₑ = sum( (-Mₗₑ*fdofs[i+2]) * BinvL[i, m] * (modofs[m+2]'*Mₗₑ) for i in 1:2, m in 1:2)

Uₗₑ = Mₗₑ*(sum(fdofs[i]*A_inL[i]-fdofs[i+2]*(αL[i,1]+αL[i,2]) for i in 1:2))

## RIGHT PORT: incoming waves 5, 6 and outgoing (transmitted) waves 7, 8
println("Building right port...")

if formulation == "SP"
    aᵣᵢ((u,p,s),(v,q,t)) = ∫( F(v,q)'⋅ (TΓ ⋅ F(u,p)) )dΓᵣ
else
    aᵣᵢ((u,p),(v,q)) = ∫( F(v,q)'⋅ (TΓ ⋅ F(u,p)) )dΓᵣ
end
Mᵣᵢ = assemble_matrix(aᵣᵢ,X,Y)

B_R = zeros(Complex{Float64}, 2, 2)
A_R = zeros(Complex{Float64}, 2, 2)
αR = zeros(Complex{Float64}, 2, 2)
for i = 1:2
    for m = 1:2
        B_R[i,m] = modofs[i+6]'*Mᵣᵢ*modofs[m+6]
        A_R[i,m] = modofs[i+6]'*Mᵣᵢ*modofs[m+4]
    end
end
BinvR = inv(B_R)

for i = 1:2
    for m = 1:2
        αR[i,m] = BinvR[i,m]*(A_R[m,1]*A_inR[1] + A_R[m,2]*A_inR[2])
    end
end

PEᵣᵢ = sum((-Mᵣᵢ*fdofs[i+6]) * BinvR[i, m] * (modofs[m+6]'*Mᵣᵢ) for i in 1:2, m in 1:2)

Uᵣᵢ = Mᵣᵢ*(sum(fdofs[i+4]*A_inR[i]-fdofs[i+6]*(αR[i,1]+αR[i,2]) for i in 1:2))
