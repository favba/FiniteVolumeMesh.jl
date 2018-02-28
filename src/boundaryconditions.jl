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