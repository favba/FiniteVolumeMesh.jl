
@inline triangle_area(n1::AbstractNode,n2::AbstractNode,n3::AbstractNode) =  norm((n2-n1) × (n3-n1))/2

@inline triangle_center(n1::AbstractNode,n2::AbstractNode,n3::AbstractNode) = (n1+n2+n3)/3

@inline normal_to_2Dline(n1::Node2D,n2::Node2D) = Node2D(n2.y - n1.y, n1.x - n2.x)
@inline normal_to_2Dline(I::NTuple{2,Integer},nodes::Vector{<:Node2D}) = normal_to_2Dline(nodes[I[1]], nodes[I[2]])

@inline line_center(n1::AbstractNode,n2::AbstractNode) = (n1+n2)/2