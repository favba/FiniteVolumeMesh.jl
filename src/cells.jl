
struct TriangleCell <: Abstract2DCell
  node1::UInt #Index of the first node
  node2::UInt #Index of the second node
  node3::UInt #Index of the third node
end

vtkcelltype(a::TriangleCell) = 5
nodetype(::Type{TriangleCell}) = Node2D{Float64}

function aspect_ratio(a::TriangleCell,m::AbstractMesh)
  l1 = distance(m.nodes[a.node1], m.nodes[a.node2])
  l2 = distance(m.nodes[a.node2], m.nodes[a.node3])
  l3 = distance(m.nodes[a.node3], m.nodes[a.node1])
  return min(l1,l2,l3)/max(l1,l2,l3)
end

volume(a::TriangleCell,n1::AbstractNode,n2::AbstractNode,n3::AbstractNode) = triangle_area(n1,n2,n3)
volume(a::TriangleCell,nodes::Vector{<:Node2D}) = volume(a, nodes[a.node1], nodes[a.node2], nodes[a.node3])
