abstract type BoundaryCondition end

struct Neumman{Q} <: BoundaryCondition end
flux(a::Neumman{Q}) where Q = Q

struct Dirichlet{Val} <: BoundaryCondition end
val(a::Dirichlet{Val}) where Val = Val

function flux_dirichlet(a::Dirichlet,bfc,bfn,k,bccenter,Tc)
  Tf = val(a)
  AfoLac = (bfn⋅bfn)/(bfn⋅(bfc - bccenter)) 
  qf = k*(Tf-Tc)*AfoLac
  return qf
end

struct FieldAtBoundary{T,BcondType,VecType}
  Tc::Array{T}
  bf2c::Array{BFace2Cell}
  bt::Array{UInt}
  bcond::BcondType
  bfc::Array{VecType}
  bccenter::Array{VecType}
  bfn::Array{VecType}
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
      b = Dirichlet{el[2]}()
    elseif el[1] == :grad
      b = Neumman{el[2]}()
    end
    bcond = (bcond...,b)
  end
  return bcond
end