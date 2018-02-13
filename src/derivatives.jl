function gradient!(out::Array{<:AbstractVec,1}, f_at_face, f_at_bface, bf2c,bfn,bcv,f2c,fn,cv,floop::Face2CellLoop{NBF,NF,VecType}) where {NBF,NF,VecType}
  fill!(out,zero(eltype(out)))
  @inbounds for i=1:NBF
    j = bf2c[i][1]
    flux = bfn[i] * f_at_bface[i] / bcv[i]
    out[j] += flux
  end  

  @inbounds for i=1:NF
    el = f2c[i]
    j1 = el[1]
    j2 = el[2]
    flux = fn[i] * f_at_face[i] 
    v = cv[i]
    out[j1] += flux/v[1]

    out[j2] -= flux/v[2]
  end  
end

@inline gradient!(out, ff, fbf, floop) = gradient!(out, ff, fbf, floop.bf2c,floop.bfn,floop.bcv,floop.f2c,floop.fn,floop.cv,floop)

function div!(out::Array{<:Number,1}, f_at_face, f_at_bface, bf2c, bfn, bcv, f2c, fn, cv, floops::Face2CellLoop{NBF,NF,VecT}) where {NBF,NF,VecT}
  fill!(out,zero(eltype(out)))
  @inbounds for i=1:NBF
    el = bf2c[i]
    j = el[1]
    n = bfn[i]
    flux = (f_at_bface[i] ⋅ n) * (1/bcv[i])
    out[j] += flux
  end  

  @inbounds for i=1:NF
    el = f2c[i]
    j1 = el[1]
    j2 = el[2]
    n = fn[i]
    flux = f_at_face[i] ⋅ n
    v = cv[i]
    out[j1] += flux * (1/v[1])
    out[j2] -= flux * (1/v[2])
  end  
end

@inline div!(out, ff, fbf, floop) = div!(out, ff, fbf, floop.bf2c,floop.bfn,floop.bcv,floop.f2c,floop.fn,floop.cv,floop)