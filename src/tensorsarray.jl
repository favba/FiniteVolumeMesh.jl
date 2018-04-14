abstract type AbstractTenArray{T<:AbstractTen,L,N} <: AbstractArray{T,N} end

@inline Base.size(u::AbstractTenArray) =
    size(u.xx)

@inline Base.length(::Type{<:AbstractTenArray{T,L,N}}) where {T,L,N} =
    L

@inline Base.length(u::AbstractTenArray) =
    length(typeof(u))

@inline Base.linearindices(::Type{<:AbstractTenArray{T,L,N}}) where {T,L,N} =
    Base.OneTo(L)

@inline Base.linearindices(u::AbstractTenArray) =
    linearindices(typeof(u))

Base.IndexStyle(u::AbstractTenArray) =
    IndexLinear()

struct Ten2DArray{T,L,N} <: AbstractTenArray{Ten2D{T},L,N}
    xx::Array{T,N}
    yx::Array{T,N}
    xy::Array{T,N}
    yy::Array{T,N}

    function Ten2DArray{T,L,N}(xx,yx,xy,yy) where {T,L,N}
        ((length(xx) == L && length(xy) == L && length(yx) == L && length(yy) == L) && 
            (size(xx) == size(yx) == size(xy) == size(yy))) || throw(ArgumentError("Arrays must have the same dimensions"))
        return new{T,L,N}(xx,yx,xy,yy)
    end

end

Ten2DArray(xx,yx,xy,yy) =
    Ten2DArray{eltype(xx),length(xx),ndims(xx)}(xx,yx,xy,yy)

@inline function Base.getindex(u::Ten2DArray{T,L,N},I) where {T,L,N}
    xx = u.xx
    yx = u.yx
    xy = u.xy
    yy = u.yy
    @boundscheck checkbounds(xx,I)
    @inbounds begin
        xxv = xx[I]
        yxv = yx[I]
        xyv = xy[I]
        yyv = yy[I]
    end
    return Ten2D{T}(xxv,yxv,xyv,yyv)
end

@inline function Base.setindex!(u::Ten2DArray{T,L,N},v::Ten2D,I) where {T,L,N}
    xx = u.xx
    yx = u.yx
    xy = u.xy
    yy = u.yy
    @boundscheck checkbounds(xx,I)
    @inbounds begin
        xx[I] = v.xx
        yx[I] = v.yx
        xy[I] = v.xy
        yy[I] = v.yy
    end
    return u
end
