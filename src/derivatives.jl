function gradient!(out::Vector, f_at_face, f_at_bface, bf2c,bfn,bcv,f2c,fn,cv,floop::Face2CellLoop)
    fill!(out,zero(eltype(out)))
    NBF = nbfaces(floop)
    @inbounds for i=1:NBF
        j = bf2c[i][1]
        flux = bfn[i] * f_at_bface[i] * bcv[i]
        out[j] += flux
    end  

    NF = nfaces(floop)
    @inbounds for i=1:NF
        j1,j2 = f2c[i]
        flux = fn[i] * f_at_face[i] 
        v = cv[i]
        out[j1] += flux*v[1]
        out[j2] -= flux*v[2]
    end  
end

@inline gradient!(out, ff, fbf, floop) = gradient!(out, ff, fbf, floop.bf2c,floop.bfn,floop.bcv,floop.f2c,floop.fn,floop.cv,floop)

@inline function gradient(ff,fbf,mesh)
    out = zeros(eltype(mesh.nodes),length(mesh.cells))
    gradient!(out,ff,fbf,mesh.f2cloops)
    return out
end

function div!(out::Vector, f_at_face, f_at_bface, bf2c, bfn, bcv, f2c, fn, cv, floops::Face2CellLoop)
    fill!(out,zero(eltype(out)))
    NBF = nbfaces(floops)
    @inbounds for i=1:NBF
        el = bf2c[i]
        j = el[1]
        n = bfn[i]
        flux = (f_at_bface[i] ⋅ n) * bcv[i]
        out[j] += flux
    end  

    NF = nfaces(floops)
    @inbounds for i=1:NF
        el = f2c[i]
        j1 = el[1]
        j2 = el[2]
        n = fn[i]
        flux = f_at_face[i] ⋅ n
        v = cv[i]
        out[j1] += flux * v[1]
        out[j2] -= flux * v[2]
    end  
end

@inline div!(out, ff, fbf, floop) = div!(out, ff, fbf, floop.bf2c,floop.bfn,floop.bcv,floop.f2c,floop.fn,floop.cv,floop)
