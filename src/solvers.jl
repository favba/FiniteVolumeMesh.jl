abstract type AbstractProblem end

struct CellProblemAdvecTemp{MeshType,TbfType,TypeK,VType,BcondType,LaplacianStruct,AdvectionStruct} <: AbstractProblem
  Tc::Array{Float64}
  k::TypeK
  ρC::Float64
  u::VType
  bcond::BcondType
  mesh::MeshType
  laplacian!::LaplacianStruct
  advection!::AdvectionStruct
  ∇Tc::Array{Vec2D{Float64}}
  Tbf::TbfType
end

function CellProblemAdvecTemp(Tc,mesh,d)
  #mesh = HomogeneousMesh(d)
  bcond = boundary_conditions(d)

  k = ConstVec{d[:conductivity]}()
  ρC = d[Symbol("rho*C")]
  uvec = ConstVec{d[:velocity]}()
  ∇Tc = zeros(Vec2D{Float64},length(Tc))
  Tbf = FieldAtBoundary(Tc,mesh,bcond)

  advection = UpWind2ndOrder()
  #laplacian = CorrectedGauss(Tc,mesh)
  laplacian = get_laplacian_method(Tc,mesh,d)
  types = typeof.((mesh,Tbf,k,uvec,bcond,laplacian,advection))
  return CellProblemAdvecTemp{types...}(Tc,k,ρC,uvec,bcond,mesh,laplacian,advection,∇Tc,Tbf)
end

function calculate_rhs!(rhs,p::CellProblemAdvecTemp)
  fill!(rhs,zero(eltype(rhs)))
  p.laplacian!(rhs,p)
  if ((typeof(p.laplacian!) == MethodB) && (typeof(p.advection!) == UpWind2ndOrder))
    gradient!(p.∇Tc,FaceSimpleInterpolation(p.Tc,p.mesh),p.Tbf,p.mesh.f2cloops)
  end
  p.advection!(rhs,p)
end