abstract type BoundaryCondition{T} end
Base.eltype(::Type{<:BoundaryCondition{T}}) where {T} = T
Base.eltype(a::BoundaryCondition) = eltype(typeof(a))

struct Neumman{T} <: BoundaryCondition{T}
    x::T
end
@inline flux(a::A) where {A<:Neumman} = a.x

struct Dirichlet{T} <: BoundaryCondition{T}
    x::T
end
@inline val(a::A) where A<:Dirichlet = a.x

function flux_dirichlet(a::Dirichlet,bfc,bfn,k,bccenter,Tc)
    Tf = val(a)
    AfoLac = (bfn⋅bfn)/(bfn⋅(bfc - bccenter)) 
    qf = k*(Tf-Tc)*AfoLac
    return qf
end

struct FieldAtBoundary{T,BcondType,VecType}
    Tc::Vector{T}
    bf2c::Vector{BFace2Cell}
    bt::Vector{UInt}
    bcond::BcondType
    bfc::Vector{VecType}
    bccenter::Vector{VecType}
    bfn::Vector{VecType}
end

function FieldAtBoundary(Tc,mesh,bcond)
    FieldAtBoundary{eltype(Tc),typeof(bcond),eltype(mesh.nodes)}(Tc,mesh.f2cloops.bf2c,mesh.f2cloops.bft,bcond,mesh.f2cloops.bfc, mesh.f2cloops.bccenter,mesh.f2cloops.bfn)
end

@inline Base.@propagate_inbounds function Base.getindex(fab::FieldAtBoundary,i)
    j = fab.bf2c[i][1]
    a = fab.bcond[fab.bt[i]] 
    return typeof(a) <: Dirichlet ? val(a) : fab.Tc[j] + ((fab.bfc[i] - fab.bccenter[i])⋅fab.bfn[i])*flux(a)/norm(fab.bfn[i])
end

function boundary_conditions(d::Dict)
    btypes = d[:BCs]
    bcond = ()
    for el in btypes
        if el[1] == :value
            b = Dirichlet(el[2])
        elseif el[1] == :grad
            b = Neumman(el[2])
        end
        bcond = (bcond...,b)
    end
    return bcond
end

function uboundary_conditions(d::Dict)
    btypes = d[:uBCs]
    bcond = ()
    for el in btypes
        if el[1] == :value
            b = Dirichlet(el[2])
        elseif el[1] == :grad
            b = Neumman(el[2])
        end
        bcond = (bcond...,b)
    end
    return bcond
end

function grad_type(bcond)
    b = eltype(bcond[1])
    if b <: Number
        return Vec2D{b}
    elseif b <: Vec2D
        return Ten2D{Float64}
    end
end
