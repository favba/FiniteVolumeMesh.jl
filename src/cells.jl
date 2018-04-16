
nodetype(a::Type{<:AbstractCell2Node}) = Vec{Float64}
@inline nodetype(a::A) where {A<:AbstractCell2Node} = nodetype(A)
nnodes(a::Type{<:AbstractCell2Node{NN,NF}}) where {NN,NF} = NN
@inline nnodes(a::A) where {A<:AbstractCell2Node} = nnodes(A)
nfaces(a::Type{<:AbstractCell2Node{NN,NF}}) where {NN,NF} = NF
@inline nfaces(a::A) where {A<:AbstractCell2Node} = nfaces(A)

nodetype(a::Type{<:AbstractCell2Node2D}) = Vec2D{Float64}
@inline nodetype(a::A) where {A<:AbstractCell2Node2D} = nodetype(A)

struct TriangleCell <: AbstractCell2Node2D{3,3}
    ind::Tuple{UInt,UInt,UInt}
end

function Base.show(io::IO,t::TriangleCell)
    compact = get(io, :compact, false)

    if !compact
        print(io,"Triangle")
    end
    print(io,"(",t[1],",",t[2],",",t[3],")") 
end

vtkcelltype(a::TriangleCell) = 5

function aspect_ratio(a::TriangleCell,m::AbstractMesh)
    @inbounds begin
        l1 = distance(m.nodes[a[1]], m.nodes[a[2]])
        l2 = distance(m.nodes[a[2]], m.nodes[a[3]])
        l3 = distance(m.nodes[a[3]], m.nodes[a[1]])
    end
    return min(l1,l2,l3)/max(l1,l2,l3)
end

cell_faces(a::TriangleCell) = Face2Node{2}((a[1], a[2])), Face2Node{2}((a[2], a[3])), Face2Node{2}((a[3], a[1]))

@inline volume(a::TriangleCell,n1::AbstractVec,n2::AbstractVec,n3::AbstractVec) = triangle_area(n1,n2,n3)
@inline function volume(a::TriangleCell,nodes::Vector{<:Vec2D}) ; @inbounds r = volume(a, nodes[a[1]], nodes[a[2]], nodes[a[3]]); return r;end

@inline perimeter(a::TriangleCell,n1::AbstractVec,n2::AbstractVec,n3::AbstractVec) = triangle_perimeter(n1,n2,n3)
@inline function perimeter(a::TriangleCell,nodes::Vector{<:Vec2D}) ; @inbounds r = perimeter(a, nodes[a[1]], nodes[a[2]], nodes[a[3]]); return r;end

@inline mean_dx(a::TriangleCell,n1::AbstractVec,n2::AbstractVec,n3::AbstractVec) = 4*volume(a,n1,n2,n3)/perimeter(a,n1,n2,n3)
@inline function mean_dx(a::TriangleCell,nodes::Vector{<:Vec2D}) ; @inbounds r = mean_dx(a, nodes[a[1]], nodes[a[2]], nodes[a[3]]); return r; end

@inline center(a::TriangleCell,n1::AbstractVec,n2::AbstractVec,n3::AbstractVec) = triangle_center(n1,n2,n3)
@inline function center(a::TriangleCell,nodes::Vector{<:Vec2D}) ; @inbounds r = center(a, nodes[a[1]], nodes[a[2]], nodes[a[3]]); return r;end

abstract type AbstractCell2Face{T,N} <: NTupleWrap{N,Int} end

#=
Idea: make lists of (volumec1,volume,c2), (centerc1,centerc2), (cellpropertyc1,cellpropertyc2)
=#
