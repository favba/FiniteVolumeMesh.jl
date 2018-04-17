struct CellField{T,VecT<:AbstractArray{T,1},BcondType} <: AbstractArray{T,1}
    f::VecT
    bcond::BcondType
end

CellField(u,bcond) =
    CellField{eltype(u),typeof(u),typeof(bcond)}(u,bcond)

function CellField(bcond,m::AbstractMesh)
    nc = length(m.cells)
    t = eltype(bcond[1])
    if t <: AbstractVec
        et = eltype(t)
        if t <: Vec2D
            u = Vec2DArray{et,1}(zeros(nc), zeros(nc))
        end
    elseif t <: Number
        u = zeros(t,nc)
    end
    CellField(u,bcond)
end

Base.size(u::CellField) =
    size(u.f)
    
Base.length(u::CellField) =
    length(u.f)

Base.IndexStyle(::Type{<:CellField}) =
    IndexLinear()

@inline function Base.getindex(u::CellField,I)
    x = u.f
    @boundscheck checkbounds(x,I)
    @inbounds xv = x[I]
    return xv
end
    
@inline function Base.setindex!(u::CellField,v,I)
    x = u.f
    @boundscheck checkbounds(x,I)
    @inbounds x[I] = v
    return u
end
    