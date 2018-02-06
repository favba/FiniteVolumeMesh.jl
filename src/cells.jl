
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

issame(a,b) = a[1] == b[2] && a[2] == b[1]

cell_faces(a::TriangleCell) = (a.node1, a.node2), (a.node2, a.node3), (a.node3, a.node1)

#function get_cell_faces(a::TriangleCell)
#  a1 = a.node1
#  a2 = a.node2
#  a3 = a.node3
#  r1 = min(a1, a2, a3)
#  r3 = max(a1, a2, a3)
#  r2 = if r1 == a1
#    if r3 == a3
#      a2
#    else
#      a3
#    end
#  elseif r1 == a2
#    if r3 == a1
#      a3
#    else
#      a1
#    end
#  elseif r1 == a3
#    if r3 == a1
#      a2
#    else
#      a1
#    end
#  end
#  return (r1,r2), (r2,r3), (r3,r1)
#end