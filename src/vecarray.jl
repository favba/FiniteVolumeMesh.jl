abstract type AbstractVecArray{T<:AbstractVec,N} <: AbstractArray{T,N} end

@inline Base.size(u::AbstractVecArray) =
    size(u.x)

@inline Base.length(u::AbstractVecArray) =
    length(u.x)

@inline Base.linearindices(u::AbstractVecArray) =
    linearindices(u.x)

Base.IndexStyle(::Type{<:AbstractVecArray}) =
    IndexLinear()

struct VecArray{T<:Number,N} <: AbstractVecArray{Vec{T},N}
    x::Array{T,N}
    y::Array{T,N}
    z::Array{T,N}

    function VecArray{T,N}(x,y,z) where {T,N}
        (size(x) == size(y) == size(z)) || throw(ArgumentError("Arrays must have the same dimensions"))
        return new{T,N}(x,y,z)
    end

end

VecArray(x,y,z) =
    VecArray{eltype(x),ndims(x)}(x,y,z)

@inline function Base.getindex(u::VecArray{T,N},I) where {T,N}
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

@inline function Base.setindex!(u::VecArray{T,N},v::Vec,I) where {T,N}
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

struct Vec2DArray{T<:Number,N} <: AbstractVecArray{Vec2D{T},N}
    x::Array{T,N}
    y::Array{T,N}

    function Vec2DArray{T,N}(x,y) where {T,N}
        (size(x) == size(y)) || throw(ArgumentError("Arrays must have the same dimensions"))
        return new{T,N}(x,y)
    end

end

Vec2DArray(x::AbstractArray,y) =
    Vec2DArray{eltype(x),ndims(x)}(x,y)

Vec2DArray{T}(L::Integer) where T = Vec2DArray{T,1}(zeros(T,L),zeros(T,L))

Base.similar(a::Vec2DArray{T}) where {T} = Vec2DArray{T}(length(a))

@inline function Base.getindex(u::Vec2DArray{T,N},I) where {T,N}
    x = u.x
    y = u.y
    @boundscheck checkbounds(x,I)
    @inbounds begin
        xv = x[I]
        yv = y[I]
    end
    return Vec2D{T}(xv,yv)
end

@inline function Base.setindex!(u::Vec2DArray{T,N},v::Vec2D,I) where {T,N}
    x = u.x
    y = u.y
    @boundscheck checkbounds(x,I)
    @inbounds begin
        x[I] = v.x
        y[I] = v.y
    end
    return u
end

struct Vec1DArray{T<:Number,N} <: AbstractVecArray{Vec1D{T},N}
    x::Array{T,N}

    function Vec1DArray{T,N}(x) where {T,N}
        return new{T,N}(x)
    end

end

Vec1DArray(x::AbstractArray) =
    Vec1DArray{eltype(x),ndims(x)}(x)

@inline function Base.getindex(u::Vec1DArray{T,N},I) where {T,N}
    x = u.x
    @boundscheck checkbounds(x,I)
    @inbounds begin
        xv = x[I]
    end
    return Vec1D{T}(xv)
end

@inline function Base.setindex!(u::Vec1DArray{T,N},v::Vec1D,I) where {T,N}
    x = u.x
    @boundscheck checkbounds(x,I)
    @inbounds begin
        x[I] = v.x
    end
    return u
end
