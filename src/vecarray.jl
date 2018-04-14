abstract type AbstractVecArray{T<:AbstractVec,L,N} <: AbstractArray{T,N} end

@inline Base.size(u::AbstractVecArray) =
    size(u.x)

@inline Base.length(::Type{<:AbstractVecArray{T,L,N}}) where {T,L,N} =
    L

@inline Base.length(u::AbstractVecArray) =
    length(typeof(u))

@inline Base.linearindices(::Type{<:AbstractVecArray{T,L,N}}) where {T,L,N} =
    Base.OneTo(L)

@inline Base.linearindices(u::AbstractVecArray) =
    linearindices(typeof(u))

Base.IndexStyle(u::AbstractVecArray) =
    IndexLinear()

struct VecArray{T<:Number,L,N} <: AbstractVecArray{Vec{T},L,N}
    x::Array{T,N}
    y::Array{T,N}
    z::Array{T,N}

    function VecArray{T,L,N}(x,y,z) where {T,L,N}
        ((length(x) == L && length(y) == L && length(z) == L) && 
            (size(x) == size(y) == size(z))) || throw(ArgumentError("Arrays must have the same dimensions"))
        return new{T,L,N}(x,y,z)
    end

end

VecArray(x,y,z) =
    VecArray{eltype(x),length(x),ndims(x)}(x,y,z)

@inline function Base.getindex(u::VecArray{T,L,N},I) where {T,L,N}
    x = u.x
    y = u.y
    z = u.z
    @boundscheck checkbounds(x,I)
    @inbounds begin
        xv = x[I]
        yv = y[I]
        zv = z[I]
    end
    return Vec{T}(xv,yv,zv)
end

@inline function Base.setindex!(u::VecArray{T,L,N},v::Vec,I) where {T,L,N}
    x = u.x
    y = u.y
    z = u.z
    @boundscheck checkbounds(x,I)
    @inbounds begin
        x[I] = v.x
        y[I] = v.y
        z[I] = v.z
    end
    return u
end

struct Vec2DArray{T<:Number,L,N} <: AbstractVecArray{Vec2D{T},L,N}
    x::Array{T,N}
    y::Array{T,N}

    function Vec2DArray{T,L,N}(x,y) where {T,L,N}
        ((length(x) == L && length(y) == L) && 
            (size(x) == size(y))) || throw(ArgumentError("Arrays must have the same dimensions"))
        return new{T,L,N}(x,y)
    end

end

Vec2DArray(x,y) =
    Vec2DArray{eltype(x),length(x),ndims(x)}(x,y)

@inline function Base.getindex(u::Vec2DArray{T,L,N},I) where {T,L,N}
    x = u.x
    y = u.y
    @boundscheck checkbounds(x,I)
    @inbounds begin
        xv = x[I]
        yv = y[I]
    end
    return Vec2D{T}(xv,yv)
end

@inline function Base.setindex!(u::Vec2DArray{T,L,N},v::Vec2D,I) where {T,L,N}
    x = u.x
    y = u.y
    @boundscheck checkbounds(x,I)
    @inbounds begin
        x[I] = v.x
        y[I] = v.y
    end
    return u
end

struct Vec1DArray{T<:Number,L,N} <: AbstractVecArray{Vec1D{T},L,N}
    x::Array{T,N}

    function Vec1DArray{T,L,N}(x) where {T,L,N}
        length(x) == L || throw(ArgumentError("Wrong lentght"))
        return new{T,L,N}(x)
    end

end

Vec1DArray(x) =
    Vec1DArray{eltype(x),length(x),ndims(x)}(x)

@inline function Base.getindex(u::Vec1DArray{T,L,N},I) where {T,L,N}
    x = u.x
    @boundscheck checkbounds(x,I)
    @inbounds begin
        xv = x[I]
    end
    return Vec1D{T}(xv)
end

@inline function Base.setindex!(u::Vec1DArray{T,L,N},v::Vec1D,I) where {T,L,N}
    x = u.x
    @boundscheck checkbounds(x,I)
    @inbounds begin
        x[I] = v.x
    end
    return u
end
