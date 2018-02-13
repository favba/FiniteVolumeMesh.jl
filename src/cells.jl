
@inline nodetype(a::Union{<:AbstractCell2Node,Type{<:AbstractCell2Node}}) = Vec{Float64}
@inline nnodes(a::Union{<:AbstractCell2Node{NN,NF},Type{<:AbstractCell2Node{NN,NF}}}) where {NN,NF} = NN
@inline nfaces(a::Union{<:AbstractCell2Node{NN,NF},Type{<:AbstractCell2Node{NN,NF}}}) where {NN,NF} = NF

@inline nodetype(a::Union{<:AbstractCell2Node2D,<:AbstractCell2Node2D}) = Vec2D{Float64}

struct TriangleCell <: AbstractCell2Node2D{3,3}
  ind::Tuple{UInt,UInt,UInt}
end

vtkcelltype(a::TriangleCell) = 5

function aspect_ratio(a::TriangleCell,m::AbstractMesh)
  l1 = distance(m.nodes[a[1]], m.nodes[a[2]])
  l2 = distance(m.nodes[a[2]], m.nodes[a[3]])
  l3 = distance(m.nodes[a[3]], m.nodes[a[1]])
  return min(l1,l2,l3)/max(l1,l2,l3)
end

cell_faces(a::TriangleCell) = Face2Node{2}((a[1], a[2])), Face2Node{2}((a[2], a[3])), Face2Node{2}((a[3], a[1]))
volume(a::TriangleCell,n1::AbstractVec,n2::AbstractVec,n3::AbstractVec) = triangle_area(n1,n2,n3)
volume(a::TriangleCell,nodes::Vector{<:Vec2D}) = volume(a, nodes[a[1]], nodes[a[2]], nodes[a[3]])

center(a::TriangleCell,n1::AbstractVec,n2::AbstractVec,n3::AbstractVec) = triangle_center(n1,n2,n3)
center(a::TriangleCell,nodes::Vector{<:Vec2D}) = center(a, nodes[a[1]], nodes[a[2]], nodes[a[3]])

abstract type AbstractCell2Face{T,N} <: NTupleWrap{N,Int} end

#=
Idea: make lists of (volumec1,volume,c2), (centerc1,centerc2), (cellpropertyc1,cellpropertyc2)
=#