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
@inline ifaces(a::MimeticGrad) = ibfaces(typeof(a))

function evalT(A::MimeticGrad,out::FaceVector{T},inp::FaceVector{T}) where {T}
  bf2c = A.f2cloop.bf2c
  brcf = A.brcf
  g = A.g
  bvkinv = A.bvkinv
  fill!(g,zero(eltype(g)))
  bcond = A.bcond
  bt = A.f2cloop.bft
  @inbounds for i in ibfaces(A)
    #j = bf2c[i][1]
    if typeof(bcond[bt[i]]) <: Neumman
      out[:b,i] = 0
    end
    #g[j] += inp[:b,i]*brcf[i]*bvkinv[i]
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
    out[:b,i] = g[j] ⋅ brcf[i]
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

function evalP(A::MimeticGrad,out::FaceVector{T},inp::FaceVector{T}) where {T}
  brcf = A.brcf
  bvkinv = A.bvkinv
  bcond = A.bcond
  bt = A.f2cloop.bft

  δ = zero(T)
  @inbounds for i in ibfaces(A)
    #if false && typeof(bcond[bt[i]]) <: Neumman
      #out[:b,i] = zero(T)
    #else
      out[:b,i] = inp[:b,i]/((brcf[i]⋅brcf[i])*bvkinv[i])
    #end
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