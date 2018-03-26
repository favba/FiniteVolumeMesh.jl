abstract type AbstractProblem end

struct CellProblemAdvecTemp{MeshType,TbfType,TypeK,TypeS,VType,BcondType,LaplacianStruct,AdvectionStruct} <: AbstractProblem
  Tc::Array{Float64}
  k::TypeK
  s::TypeS
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
  bcond = boundary_conditions(d)

  k = ConstVec{d[:conductivity]}()
  s = ConstVec{d[:source]}()
  ρC = d[Symbol("rho*C")]
  uvec = ConstVec{d[:velocity]}()
  ∇Tc = zeros(Vec2D{Float64},length(Tc))
  Tbf = FieldAtBoundary(Tc,mesh,bcond)

  advection = get_advection_method(d)
  laplacian = get_laplacian_method(Tc,bcond,k,mesh,d)
  types = typeof.((mesh,Tbf,k,s,uvec,bcond,laplacian,advection))
  return CellProblemAdvecTemp{types...}(Tc,k,s,ρC,uvec,bcond,mesh,laplacian,advection,∇Tc,Tbf)
end

@inline needs_gradient_calculation(p) = (typeof(p.advection!) === UpWind2ndOrder) && (typeof(p.laplacian!) === MethodB) && (typeof(p.u) <: ConstVec && p.u[1] == zero(eltype(p.u)))

function calculate_rhs!(rhs,p::CellProblemAdvecTemp)
  fill!(rhs,zero(eltype(rhs)))

  is_implicit(p.laplacian!) || p.laplacian!(rhs,p)

  needs_gradient_calculation(p) && gradient!(p.∇Tc,FaceSimpleInterpolation(p.Tc,p.mesh),p.Tbf,p.mesh.f2cloops)

  if !(typeof(p.u) <: ConstVec && p.u[1] == zero(eltype(p.u))) 
    p.advection!(rhs,p)
  end
  if !(typeof(p.s) <: ConstVec && p.s[1] == zero(eltype(p.s))) 
    add_source!(rhs,p)
  end

  is_implicit(p.laplacian!) && p.laplacian!(rhs,p)

  return nothing
end

function add_source!(rhs,p)
  s = p.s
  @inbounds @simd for i=1:length(rhs)
    rhs[i] += s[i]
  end
end