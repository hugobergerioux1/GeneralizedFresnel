# Differential operators of the 2D mixed formulation (fields independent of y).
#   Gridap coordinates: x[1] = z (normal to the interface), x[2] = x (along the interface)
#   unknowns: u = (Ez, Ex) in H(curl) (Nédélec), p = Ey in H¹ (Lagrange)

# partial derivatives of a scalar field p, and of the component i of a vector field u
∂z(p::CellField) = Operation(_z)(∇(p))
∂x(p::CellField) = Operation(_x)(∇(p))
_z(u::VectorValue) = u[1]
_x(u::VectorValue) = u[2]
∂z(u::CellField,i::Int64) = Operation(col1)(∇(u),i)
∂x(u::CellField,i::Int64) = Operation(col2)(∇(u),i)
col1(M::TensorValue,i) = M[1,i]
col2(M::TensorValue,i) = M[2,i]

v1 = VectorValue(1,0,0)
v2 = VectorValue(0,1,0)
v3 = VectorValue(0,0,1)

# 3D field E = (Ex, Ey, Ez) from (u, p)
F(u::CellField,p::CellField) = Operation(newF)(u,p)
newF(u::VectorValue,p) = VectorValue(u[2],p,u[1])

# ∇×E = (−∂z Ey, ∂z Ex − ∂x Ez, ∂x Ey)
curl2D(u::CellField,p::CellField) = -∂z(p)              *v1 +
                                    (∂z(u,2)-∂x(u,1))   *v2 +
                                    ∂x(p)               *v3

# ∇s = (∂x s, 0, ∂z s) of a scalar field s (Lagrange multiplier of the SP formulation)
function grad2D(f)
    return ∂x(f) * v1 + ∂z(f) * v3
end

# component-wise modulus and phase of a complex vector field
MAG(E::CellField) = Operation(vmag)(E)
vmag(E::VectorValue) = VectorValue(abs.(Tuple(E))...)
phase(E::CellField) = Operation(vangle)(E)
vangle(E::VectorValue) = VectorValue(angle.(Tuple(E))...)
