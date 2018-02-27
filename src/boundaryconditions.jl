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