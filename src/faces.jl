
@inline area(face::AbstractFace2Node1D,nodes::Vector{<:Vec1D}) = 1 
@inline area(face::AbstractFace2Node2D,nodes::Vector{<:Vec2D}) =  normal_to_2Dline(face.ind,nodes)

@inline function center(face::AbstractFace2Node1D,nodes::Vector{<:Vec1D}); @inbounds r = nodes[face[1]]; return r; end 
@inline center(face::AbstractFace2Node2D,nodes::Vector{<:Vec2D}) =  line_center(face.ind,nodes)

struct Face2Node{N} <: AbstractFace2Node{N}
    ind::NTuple{N,UInt}
end

function Base.show(io::IO,f::Face2Node{N}) where N
    compact = get(io, :compact, false)

    if !compact
        print(io,"$(N)DFace")
    end
    print(io,"(",join(f.ind,","),")") 
end

abstract type AbstractFace2Cell <: NTupleWrap{2,UInt} end

abstract type AbstractBFace2Cell <: NTupleWrap{1,UInt} end

struct Face2Cell <: AbstractFace2Cell
    ind::Tuple{UInt,UInt}
end

function Base.show(io::IO,f::Face2Cell)
    compact = get(io, :compact, false)

    if !compact
        print(io,"Interface of cells ",f[1]," and ",f[2])
    else
        print(io,"(",join(f.ind,","),")") 
    end
end

struct BFace2Cell <: AbstractBFace2Cell
    ind::Tuple{UInt}
end

function Base.show(io::IO,f::BFace2Cell)
    compact = get(io, :compact, false)

    if !compact
        print(io,"Boundary face of cell ",f[1])
    else
        print(io,f[1]) 
    end
end

struct Face2CellLoop{NBF,NF,VecType}
    bf2c::Vector{BFace2Cell}
    bft::Vector{UInt} # List of of boundary faces types
    bfn::Vector{VecType} # List of normals of boundary faces
    bfc::Vector{VecType} # List of centers of boundary faces
    bcv::Vector{Float64} # List of 1/volumes of boundary faces owner cells
    bccenter::Vector{VecType} # List of centers of boundary faces owner cells
    f2c::Vector{Face2Cell}
    fn::Vector{VecType} # List of normals of faces
    fc::Vector{VecType} # List of centers of faces
    cv::Vector{NTuple{2,Float64}} # List of 1/volumes of faces owner cells
    ccenter::Vector{NTuple{2,VecType}} # List of centers of faces owner cells
end

nbfaces(::Type{Face2CellLoop{NBF,NF,Vectype}}) where {NBF,NF,Vectype} = NBF
@inline nbfaces(a::A) where A<:Face2CellLoop = nbfaces(A)
nfaces(::Type{Face2CellLoop{NBF,NF,Vectype}}) where {NBF,NF,Vectype} = NF
@inline nfaces(a::A) where A<:Face2CellLoop = nfaces(A)

struct FaceSimpleInterpolation{T,NF} <: AbstractVector{T}
    Tc::Vector{T}
    f2c::Vector{Face2Cell}
end

FaceSimpleInterpolation(Tc,mesh) = FaceSimpleInterpolation{eltype(Tc),length(mesh.f2cloops.f2c)}(Tc,mesh.f2cloops.f2c)

Base.size(a::Type{FaceSimpleInterpolation{T,NF}}) where {T,NF} = (NF,)
@inline Base.size(a::A) where {A<:FaceSimpleInterpolation} = size(A)
Base.length(a::Type{FaceSimpleInterpolation{T,NF}}) where {T,NF} = NF
@inline Base.length(a::A) where {A<:FaceSimpleInterpolation} = length(A)
Base.linearindices(a::FaceSimpleInterpolation{T,NF}) where {T,NF} = Base.OneTo(length(a))
Base.IndexStyle(::Type{FaceSimpleInterpolation}) = IndexLinear()
@inline Base.@propagate_inbounds function Base.getindex(a::FaceSimpleInterpolation,I)
    j1,j2 = a.f2c[I]
    Tc = a.Tc
    return 0.5*(Tc[j1]+Tc[j2])
end

struct FaceVector{T,NF,NBF,L} <: AbstractVector{T}
    data::Vector{T}
end

Base.length(::Type{<:FaceVector{T,NF,NBF,L}}) where {T,NF,NBF,L} = L
Base.size(::Type{<:FaceVector{T,NF,NBF,L}}) where {T,NF,NBF,L} = (L,)
Base.linearindices(::Type{<:FaceVector{T,NF,NBF,L}}) where {T,NF,NBF,L} = Base.OneTo(L)
Base.IndexStyle(::Type{<:FaceVector{T,NF,NBF,L}}) where {T,NF,NBF,L} = IndexLinear()
Base.similar(a::FaceVector{T,NF,NBF,L}) where {T,NF,NBF,L} = FaceVector{T,NF,NBF,L}(zeros(T,L))

for op in (:(Base.length),:(Base.size),:(Base.linearindices))
    @eval begin
        @inline function ($op)(a::FaceVector)
            return ($op)(typeof(a)) 
        end
    end
end

@inline function Base.getindex(a::FaceVector,i::Integer)
    d = a.data
    @boundscheck checkbounds(a,i)
    @inbounds r = d[i]
    return r
end

@inline function Base.setindex!(a::FaceVector,x,i::Integer)
    d = a.data
    @boundscheck checkbounds(a,i)
    @inbounds d[i] = x
    return d
end

@inline @generated function Base.getindex(a::FaceVector{T,NF,NBF,L},s::Symbol,i1::Integer) where {T,NF,NBF,L}
    n = NF
    quote
        d = a.data
        i = i1 + $n
        @boundscheck checkbounds(a,i)
        @inbounds r = d[i]
        return r
    end
end

@inline @generated function Base.setindex!(a::FaceVector{T,NF,NBF,L},x,s::Symbol,i1::Integer) where {T,NF,NBF,L}
    n = NF
    quote
        d = a.data
        i = i1 + $n
        @boundscheck checkbounds(a,i)
        @inbounds d[i] = x
        return d
    end
end
