abstract type AbstractMatrixLike end

struct MimeticGrad{T,VecType,BcondType,NBF,NF} <: AbstractMatrixLike
  g::Vector{VecType} # cell stuff
  brcf::Vector{VecType} # cell-to-face vector
  rcf::Vector{NTuple{2,VecType}} # cell-to-face vector
  bvkinv::Vector{T} # 1/(cell_volume*K)
  vkinv::Vector{NTuple{2,T}} # 1/(cell_volume*K)
  bcond::BcondType
  f2cloop::Face2CellLoop{NBF,NF,VecType}
end

function MimeticGrad(k,bcond,m::HomogeneousMesh)
  T = eltype(k)
  vectype = vec_type(m)
  f2cloop = m.f2cloops
  g = zeros(vectype,length(m.cells))
  brcf, rcf, bvkinv, vkinv = get_stuff_needed_for_MimeticGrad(k,m)
  return MimeticGrad{T,vectype,typeof(bcond),nbfaces(m),nfaces(m)}(g,brcf,rcf,bvkinv,vkinv,bcond,f2cloop)
end

function get_stuff_needed_for_MimeticGrad(k,m::HomogeneousMesh)
  nbf = nbfaces(m)
  nf = nfaces(m)

  brcf = zeros(vec_type(m),nbf)
  rcf = Array{NTuple{2,vec_type(m)}}(nf)
  T = eltype(k) 
  bvkinv = zeros(T,nbf)
  vkinv = Array{NTuple{2,T}}(nf)

  nodes = m.nodes
  cells = m.cells
  bf2c = m.f2cloops.bf2c
  bfaces = m.bfaces
  bcv = m.f2cloops.bcv
  @inbounds for i in 1:nbf
    j = bf2c[i][1]
    ccenter = center(cells[j],nodes)
    brcf[i] = center(bfaces[i],nodes) - ccenter
    bvkinv[i] = bcv[i]/k[i]
  end
  
  f2c = m.f2cloops.f2c
  faces = m.faces
  cv = m.f2cloops.cv
  @inbounds for i in 1:nf
    j1,j2 = f2c[i]
    cc1 = center(cells[j1],nodes)
    cc2 = center(cells[j2],nodes)
    fc = center(faces[i],nodes) 
    rcf[i] = (fc - cc1, fc - cc2)
    cv1,cv2 = cv[i]
    vkinv[i] = (cv1/k[j1], cv2/k[j2])
  end

  return brcf, rcf, bvkinv, vkinv
end

ibfaces(::Type{MimeticGrad{T,VecType,NC,NBF,NF}}) where {T,VecType,NC,NBF,NF} = Base.OneTo(NBF)
@inline ibfaces(a::MimeticGrad) = ibfaces(typeof(a))

ifaces(::Type{MimeticGrad{T,VecType,NC,NBF,NF}}) where {T,VecType,NC,NBF,NF} = Base.OneTo(NF)
@inline ifaces(a::MimeticGrad) = ifaces(typeof(a))

function evalT(A::MimeticGrad,out::FaceVector{T},inp::FaceVector{T}) where {T}
  bf2c = A.f2cloop.bf2c
  brcf = A.brcf
  g = A.g
  bvkinv = A.bvkinv
  fill!(g,zero(eltype(g)))
  bcond = A.bcond
  bt = A.f2cloop.bft
  @inbounds for i in ibfaces(A)
    j = bf2c[i][1]
    if typeof(bcond[bt[i]]) <: Neumman
      inp[:b,i] = zero(T)
    else
      g[j] += inp[:b,i]*brcf[i]*bvkinv[i]
    end
  end
  
  f2c = A.f2cloop.f2c
  rcf = A.rcf
  vkinv = A.vkinv
  @inbounds for i in ifaces(A)
    j1,j2 = f2c[i]
    rc1,rc2 = rcf[i]
    v1,v2 = vkinv[i]
    g[j1] += inp[i]*rc1*v1
    g[j2] -= inp[i]*rc2*v2
  end
 
  δ = zero(T)
  @inbounds for i in ibfaces(A)
    j = bf2c[i][1]
    if typeof(bcond[bt[i]]) <: Neumman
      out[:b,i] = zero(T)
    else
      out[:b,i] = g[j] ⋅ brcf[i]
    end
    δ += inp[:b,i]*out[:b,i]
  end

  @inbounds for i in ifaces(A)
    j1,j2 = f2c[i]
    rc1,rc2 = rcf[i]
    out[i] = (g[j1] ⋅ rc1) - (g[j2] ⋅ rc2)
    δ += inp[i]*out[i]
  end

  return δ
end

function evalP(A::MimeticGrad,out::FaceVector{T},inp::FaceVector{T}) where {T}
  brcf = A.brcf
  bvkinv = A.bvkinv
  bcond = A.bcond
  bt = A.f2cloop.bft

  δ = zero(T)
  @inbounds for i in ibfaces(A)
    if typeof(bcond[bt[i]]) <: Neumman
      out[:b,i] = zero(T)
      inp[:b,i] = zero(T)
    else
      out[:b,i] = inp[:b,i]/((brcf[i]⋅brcf[i])*bvkinv[i])
    end
    δ += out[:b,i] * inp[:b,i]
  end
  
  rcf = A.rcf
  vkinv = A.vkinv
  @inbounds for i in ifaces(A)
    rc1,rc2 = rcf[i]
    v1,v2 = vkinv[i]
    out[i] = inp[i]/((rc1⋅rc1)*v1 + (rc2⋅rc2)*v2)
    δ += out[i] * inp[i]
  end
 
  return δ
end

function Base.A_mul_B!(out::FaceVector{T},A::MimeticGrad,inp::FaceVector{T}) where {T}
  bf2c = A.f2cloop.bf2c
  brcf = A.brcf
  g = A.g
  bvkinv = A.bvkinv
  fill!(g,zero(eltype(g)))

  @inbounds for i in ibfaces(A)
    j = bf2c[i][1]
    g[j] += inp[:b,i]*brcf[i]*bvkinv[i]
  end
  
  f2c = A.f2cloop.f2c
  rcf = A.rcf
  vkinv = A.vkinv
  @inbounds for i in ifaces(A)
    j1,j2 = f2c[i]
    rc1,rc2 = rcf[i]
    v1,v2 = vkinv[i]
    g[j1] += inp[i]*rc1*v1
    g[j2] -= inp[i]*rc2*v2
  end
 
  @inbounds for i in ibfaces(A)
    j = bf2c[i][1]
    out[:b,i] = g[j] ⋅ brcf[i]
  end

  @inbounds for i in ifaces(A)
    j1,j2 = f2c[i]
    rc1,rc2 = rcf[i]
    out[i] = (g[j1] ⋅ rc1) - (g[j2] ⋅ rc2)
  end

  return out
end

function Base.:*(A::AbstractMatrixLike,x::AbstractVector{T}) where {T}
  b = similar(x)
  return A_mul_B!(b,A,x)
end

struct ImplicitMethodB{MeshType,BcondType,Ktype,T} <: AbstractMatrixLike
  m::MeshType
  bcond::BcondType
  k::Ktype
  dtrhoc::T
end

function evalT(A::ImplicitMethodB,out::AbstractArray{T},inp::AbstractArray{T}) where {T}

  fill!(out,zero(T))

  floops = A.m.f2cloops
  bf2c = floops.bf2c
  bcond = A.bcond
  bt = floops.bft 
  bfc = floops.bfc 
  bfn = floops.bfn
  k = A.k
  bccenter = floops.bccenter
  bcv = floops.bcv

  NBF = nbfaces(floops)
   #@inbounds for i=1:NBF
    #j = bf2c[i][1]
    #bch = bcond[bt[i]] 
    #if typeof(bch) <: Neumman
      #rhs[j] += flux(bch)*bcv[i]
    #else
      #rhs[j] += x[:b,i]*bcv[i]
    #end
  #   j = bf2c[i][1]
  #   qf = bpart(bcond[bt[i]],i,j,bfc,bfn,k,bccenter,inp)
  #   out[j] += qf*bcv[i]
  #end  
   @inbounds for i=1:NBF
    j = bf2c[i][1]
    bch = bcond[bt[i]] 
    if typeof(bch) <: Neumman
      qf = flux(bch)
    else
      AfoLac = (bfn[i]⋅bfn[i])/(bfn[i]⋅(bfc[i] - bccenter[i])) 
      qf = -k[j]*inp[j]*AfoLac
    end
    out[j] += qf*bcv[i]
  end  

  f2c = floops.f2c 
  fn = floops.fn
  ccenter = floops.ccenter
  cv = floops.cv

  NF = nfaces(floops)
  @inbounds for i=1:NF
    el = f2c[i]
    j1 = el[1]
    j2 = el[2]
    To = inp[j1]
    Ta = inp[j2]
    n = fn[i] 
    c = ccenter[i]
    AfoLac = (n⋅n)/(n⋅(c[2] - c[1])) 
    qf = k[i]*(Ta-To)*AfoLac
  
    v = cv[i]
    out[j1] += qf*v[1]
    out[j2] -= qf*v[2]
  end  

  δ = zero(T)
  dtrhoc = A.dtrhoc
  @inbounds for i in linearindices(out)
    out[i] = inp[i] - dtrhoc*out[i]
    δ += out[i]*inp[i]
  end
  
  return δ
end

function evalP(A::ImplicitMethodB,out::AbstractArray{T},inp::AbstractArray{T}) where {T}

  fill!(out,zero(T))

  floops = A.m.f2cloops
  bf2c = floops.bf2c
  bcond = A.bcond
  bt = floops.bft 
  bfc = floops.bfc 
  bfn = floops.bfn
  k = A.k
  bccenter = floops.bccenter
  bcv = floops.bcv
  
  NBF = nbfaces(floops)
  @inbounds for i=1:NBF
    j = bf2c[i][1]
    bch = bcond[bt[i]] 
    if typeof(bch) <: Neumman
      qf = flux(bch)
    else
      AfoLac = (bfn[i]⋅bfn[i])/(bfn[i]⋅(bfc[i] - bccenter[i])) 
      qf = -k[j]*AfoLac
    end
    out[j] += qf*bcv[i]
  end  

  f2c = floops.f2c 
  fn = floops.fn
  ccenter = floops.ccenter
  cv = floops.cv

  NF = nfaces(floops)
  @inbounds for i=1:NF
    el = f2c[i]
    j1 = el[1]
    j2 = el[2]
    To = inp[j1]
    Ta = inp[j2]
    n = fn[i] 
    c = ccenter[i]
    AfoLac = (n⋅n)/(n⋅(c[2] - c[1])) 
    #qf = k[i]*(Ta-To)*AfoLac
  
    v = cv[i]
    out[j1] -= k[j1]*AfoLac*v[1]
    out[j2] -= k[j2]*AfoLac*v[2]
  end  

  δ = zero(T)
  dtrhoc = A.dtrhoc
  @inbounds for i in linearindices(out)
    out[i] = inp[i]/(1 - dtrhoc*out[i])
    δ += out[i]*inp[i]
  end
  
  return δ
end

function Base.A_mul_B!(out::AbstractArray{T},A::ImplicitMethodB,inp::AbstractArray{T}) where {T}
  floops = A.m.f2cloops
  bf2c = floops.bf2c
  bcond = A.bcond
  bt = floops.bft 
  bfc = floops.bfc 
  bfn = floops.bfn
  k = A.k
  bccenter = floops.bccenter
  bcv = floops.bcv
  
  fill!(out,zero(T))

  NBF = nbfaces(floops)
   @inbounds for i=1:NBF
     j = bf2c[i][1]
     qf = bpart(bcond[bt[i]],i,j,bfc,bfn,k,bccenter,inp)
     out[j] += qf*bcv[i]
   end  

  f2c = floops.f2c 
  fn = floops.fn
  ccenter = floops.ccenter
  cv = floops.cv

  NF = nfaces(floops)
  @inbounds for i=1:NF
    el = f2c[i]
    j1 = el[1]
    j2 = el[2]
    To = inp[j1]
    Ta = inp[j2]
    n = fn[i] 
    c = ccenter[i]
    AfoLac = (n⋅n)/(n⋅(c[2] - c[1])) 
    qf = k[i]*(Ta-To)*AfoLac
  
    v = cv[i]
    out[j1] += qf*v[1]
    out[j2] -= qf*v[2]
  end  

  dtrhoc = A.dtrhoc
  @inbounds for i in linearindices(out)
    out[i] = inp[i] - dtrhoc*out[i]
  end
  
  return nothing
end