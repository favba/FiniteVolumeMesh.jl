
function bpart(btype::A,i::Integer,j::Integer,bfc,bfn,k,bccenter,cfield) where {A<:BoundaryCondition}
  if A<:Neumman
    return flux(btype)
  elseif A<:Dirichlet
    return flux_dirichlet(btype,bfc[i],bfn[i],k[i],bccenter[i],cfield[j])
  end
end

# Method B
function laplacian_mB!(out,cfield,k,bf2c,bt,bcond,bfn,bfc,bccenter,bcv,f2c,fn,ccenter,cv,floops::Face2CellLoop)
  #fill!(out,zero(eltype(out)))
  NBF = nbfaces(floops)
  @inbounds for i=1:NBF
    j = bf2c[i][1]
    qf = bpart(bcond[bt[i]],i,j,bfc,bfn,k,bccenter,cfield)
    out[j] += qf*bcv[i]
  end  

  NF = nfaces(floops)
  @inbounds for i=1:NF
    el = f2c[i]
    j1 = el[1]
    j2 = el[2]
    To = cfield[j1]
    Ta = cfield[j2]
    n = fn[i] 
    c = ccenter[i]
    AfoLac = (n⋅n)/(n⋅(c[2] - c[1])) 
    qf = k[i]*(Ta-To)*AfoLac
  
    v = cv[i]
    out[j1] += qf*v[1]
    out[j2] -= qf*v[2]
  end  
end

@inline laplacian_mB!(out,cfield,k,bcond,floops) = laplacian_mB!(out,cfield,k,
  floops.bf2c,floops.bft,bcond,floops.bfn,floops.bfc,floops.bccenter,floops.bcv,
  floops.f2c,floops.fn,floops.ccenter,floops.cv,
floops)

function fillbvalues!(Tbf,bcond,Tc,bf2c,bt,bfc,bccenter,bfn,floops::Face2CellLoop{NBF,NF,VT}) where {NBF,NF,VT}
  @inbounds for i=1:NBF
    j = bf2c[i][1]
    a = bcond[bt[i]] 
    Tbf[i] = typeof(a) <: Dirichlet ? val(a) : Tc[j] + ((bfc[i] - bccenter[i])⋅bfn[i])*flux(a)/norm(bfn[i])
  end   
end

# Gauss Method
function laplacian_Gauss!(out, ∇Tc, Tc, Tf, Tbf, k,bf2c,bt,bcond,bfn,bfc,bccenter,bcv,f2c,fn,fc,ccenter,cv,floops::Face2CellLoop{NBF,NF,VT}) where {NBF,NF,VT}
  #fill!(out,zero(eltype(out)))

  #fillbvalues!(Tbf,bcond,Tc,bf2c,bt,bfc,bccenter,bfn,floops)

  #Qf at bfaces
  @inbounds for i=1:NBF
    j = bf2c[i][1]
    qf = bpart(bcond[bt[i]],i,j,bfc,bfn,k,bccenter,Tc)
    out[j] += qf*bcv[i]
  end  

  #first Tf
  if !(typeof(Tf)<:FaceSimpleInterpolation)
    @inbounds for i=1:NF
      j1,j2 = f2c[i].ind
      Tf[i] = 0.5*(Tc[j1] + Tc[j2])
    end
  end

  gradient!(∇Tc, Tf, Tbf,floops)
  
  @inbounds for i=1:NF
    j1, j2 = f2c[i]
    n = fn[i]
    v1, v2 = cv[i]
    qf = 0.5 * k[i] * (n ⋅ (∇Tc[j1] + ∇Tc[j2]))
    out[j1] += qf*v1
    out[j2] -= qf*v2
  end  

end

@inline laplacian_Gauss!(out, ∇Tc, Tc, Tf, Tbf,k,bcond,floops) = laplacian_Gauss!(out, ∇Tc, Tc, Tf, Tbf,k,
  floops.bf2c,floops.bft,bcond,floops.bfn,floops.bfc,floops.bccenter,floops.bcv,
  floops.f2c,floops.fn,floops.fc,floops.ccenter,floops.cv,
floops)

# Corrected Gauss Method
function laplacian_CGauss!(out, ∇Tc, Tc, Tf::AbstractVector, Tbf, k,bf2c,bt,bcond,bfn,bfc,bccenter,bcv,f2c,fn,fc,ccenter,cv,floops::Face2CellLoop{NBF,NF,VT}) where {NBF,NF,VT}
  #fill!(out,zero(eltype(out)))

  #fillbvalues!(Tbf,bcond,Tc,bf2c,bt,bfc,bccenter,bfn,floops)

  #Qf at bfaces
  @inbounds for i=1:NBF
    j = bf2c[i][1]
    qf = bpart(bcond[bt[i]],i,j,bfc,bfn,k,bccenter,Tc)
    out[j] += qf*bcv[i]
  end  

  #first Tf
  if !(typeof(Tf)<:FaceSimpleInterpolation)
    @inbounds for i=1:NF
      j1,j2 = f2c[i].ind
      Tf[i] = 0.5*(Tc[j1] + Tc[j2])
    end
  end
  gradient!(∇Tc, Tf, Tbf,floops)
  
  @inbounds for i=1:NF
    j1, j2 = f2c[i]
    n = fn[i]
    v1, v2 = cv[i]
    ∇T̂f = 0.5*(∇Tc[j1] + ∇Tc[j2]) 
    Lac = ccenter[i][2] - ccenter[i][1]
    L2 = Lac⋅Lac
    qf = k[i] * ((Tc[j2] - Tc[j1])*(Lac⋅n)/L2 + ∇T̂f⋅n - (∇T̂f⋅Lac)*(Lac⋅n)/L2)
    out[j1] += qf*v1
    out[j2] -= qf*v2
  end  

end

@inline laplacian_CGauss!(out, ∇Tc, Tc, Tf, Tbf,k,bcond,floops) = laplacian_CGauss!(out, ∇Tc, Tc, Tf, Tbf,k,
  floops.bf2c,floops.bft,bcond,floops.bfn,floops.bfc,floops.bccenter,floops.bcv,
  floops.f2c,floops.fn,floops.fc,floops.ccenter,floops.cv,
floops)

abstract type AbstractLaplacian end

struct CorrectedGauss{ArrayType} <: AbstractLaplacian
  Tf::ArrayType
end

CorrectedGauss(Tc,mesh) = CorrectedGauss(FaceSimpleInterpolation(Tc,mesh))

function (l::CorrectedGauss)(rhs,p)
  laplacian_CGauss!(rhs, p.∇Tc, p.Tc, l.Tf, p.Tbf, p.k, p.bcond, p.mesh.f2cloops)
end

struct Gauss{ArrayType} <: AbstractLaplacian
  Tf::ArrayType
end

Gauss(Tc,mesh) = Gauss(FaceSimpleInterpolation(Tc,mesh))

function (l::Gauss)(rhs,p)
  laplacian_Gauss!(rhs, p.∇Tc, p.Tc, l.Tf, p.Tbf, p.k, p.bcond, p.mesh.f2cloops)
end

struct MethodB <: AbstractLaplacian
end

function (l::MethodB)(rhs,p)
  laplacian_mB!(rhs, p.Tc, p.k, p.bcond, p.mesh.f2cloops)
end

struct CellMimetic{M<:AbstractMatrixLike,P,ArrayT} <: AbstractLaplacian
  A::M
  pcg::P
  b::ArrayT
  x::ArrayT
end

function CellMimetic(k,bcond,m)
  nbf = nbfaces(m)
  nf = nfaces(m)
  nt = nbf+nf
  x = FaceVector{eltype(k),nf,nbf,nt}(zeros(eltype(k),nt))
  b = similar(x)
  pcg = PCG(x)
  A = MimeticGrad(k,bcond,m)
  return CellMimetic{typeof(A),typeof(pcg),typeof(b)}(A,pcg,b,x)
end

function (l::CellMimetic)(rhs,p)
  set_x_and_b!(l,p,p.mesh.f2cloops)
  solve!(l.pcg,l.A,l.x,l.b,true)
  laplacian_mimetic!(rhs, l.x, p.bcond, p.Tc, p.∇Tc, l.A.rcf, l.A.brcf, p.k, p.mesh.f2cloops)
end

function set_x_and_b!(l,p,f2cl::Face2CellLoop)

  nbf = nbfaces(f2cl)
  nf = nfaces(f2cl)
  Tbf = p.Tbf
  Tc = p.Tc
  b = l.b
  x = l.x
  bf2c = f2cl.bf2c

  bcond = p.bcond
  bfc = f2cl.bfc
  bfn = f2cl.bfn
  bccenter = f2cl.bccenter
  k = p.k
  bt = f2cl.bft
  bcv = f2cl.bcv
  @inbounds for i = 1:nbf
    j = bf2c[i][1]
    b[:b,i] = Tbf[i] - Tc[j]
    x[:b,i] = bpart(bcond[bt[i]],i,j,bfc,bfn,k,bccenter,Tc)
  end

  f2c = f2cl.f2c

  @inbounds for i = 1:nf
    j1,j2 = f2c[i]
    b[i] = Tc[j2] - Tc[j1]
  end

end


# The gradient calculation is not right. It seems to be off by some multiplication factor.
function laplacian_mimetic!(rhs,x,bcond,Tc,∇Tc,rcf,brcf,k,f2cl::Face2CellLoop)
  nbf = nbfaces(f2cl)
  nf = nfaces(f2cl)

  bf2c = f2cl.bf2c

  bt = f2cl.bft
  bcv = f2cl.bcv
  @inbounds for i=1:nbf
    j = bf2c[i][1]
    bch = bcond[bt[i]] 
    if typeof(bch) <: Neumman
      rhs[j] += flux(bch)*bcv[i]
      ∇Tc[j] += flux(bch)*brcf[i]/k[i]
    else
      rhs[j] += x[:b,i]*bcv[i]
      ∇Tc[j] += x[:b,i]*brcf[i]/k[i]
    end
  end  

  f2c = f2cl.f2c
  cv = f2cl.cv
  @inbounds for i=1:nf
    j1, j2 = f2c[i]
    cv1,cv2 = cv[i]
    qf = x[i]
    rhs[j1] += qf*cv1
    rhs[j2] -= qf*cv2

    rc1,rc2 = rcf[i]
    ∇Tc[j1] += qf*rc1/k[i]
    ∇Tc[j2] -= qf*rc2/k[i]
  end

end

@inline is_implicit(a::AbstractLaplacian) = false
abstract type AbstractImplicitLaplacian <: AbstractLaplacian end
@inline is_implicit(a::AbstractImplicitLaplacian) = true


struct ImplicitDiff{M<:AbstractMatrixLike,P,T<:Number} <: AbstractImplicitLaplacian
  A::M
  pcg::P
  b::Vector{T}
end

function ImplicitDiff(k,bcond,mesh,d)
  dtrhoc = d[:dt]*d[Symbol("rho*C")]
  A = ImplicitMethodB{typeof(mesh),typeof(bcond),typeof(k),typeof(dtrhoc)}(mesh,bcond,k,dtrhoc)
  pcg = PCG(zeros(length(mesh.cells)))
  b = zeros(eltype(k),length(mesh.cells))
  return ImplicitDiff{typeof(A),typeof(pcg),eltype(k)}(A,pcg,b)
end

function (l::ImplicitDiff)(rhs,p)
  set_x_and_b!(l,l.b,rhs,p)
  solve!(l.pcg,l.A,rhs,l.b,true)
  return nothing
end

function set_x_and_b!(l::ImplicitDiff,b,rhs,p)
  A_mul_B!(b,l.A,p.Tc)
  Tc = p.Tc
  dtrhoc = l.A.dtrhoc
  @inbounds @simd for i in linearindices(b)
    b[i] = Tc[i] + dtrhoc*rhs[i] - b[i]
  end
  return nothing
end