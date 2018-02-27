
function bpart(btype::A,i::Integer,j::Integer,bfc,bfn,k,bccenter,cfield) where {A<:BoundaryCondition}
  if A<:Neumman
    return flux(btype)
  elseif A<:Dirichlet
    return flux_dirichlet(btype,bfc[i],bfn[i],k[i],bccenter[i],cfield[j])
  end
end

# Method B
function laplacian_mB!(out,cfield,k,bf2c,bt,bcond,bfn,bfc,bccenter,bcv,f2c,fn,ccenter,cv,floops::Face2CellLoop{NBF,NF,VT}) where {NBF,NF,VT}
  fill!(out,zero(eltype(out)))
  @inbounds for i=1:NBF
    j = bf2c[i][1]
    qf = bpart(bcond[bt[i]],i,j,bfc,bfn,k,bccenter,cfield)
    out[j] += qf*bcv[i]
  end  

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
  fill!(out,zero(eltype(out)))

  fillbvalues!(Tbf,bcond,Tc,bf2c,bt,bfc,bccenter,bfn,floops)

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
  fill!(out,zero(eltype(out)))

  fillbvalues!(Tbf,bcond,Tc,bf2c,bt,bfc,bccenter,bfn,floops)

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

