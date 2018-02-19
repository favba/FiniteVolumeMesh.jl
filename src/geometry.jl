
@inline triangle_area(n1::AbstractVec,n2::AbstractVec,n3::AbstractVec) =  0.5*norm((n2-n1) × (n3-n1))

@inline triangle_center(n1::AbstractVec,n2::AbstractVec,n3::AbstractVec) = (n1+n2+n3)/3

@inline triangle_perimeter(n1::AbstractVec,n2::AbstractVec,n3::AbstractVec) = norm(n2-n1)+norm(n3-n2)+norm(n1-n3)

@inline normal_to_2Dline(n1::Vec2D,n2::Vec2D) = Vec2D(n2.y - n1.y, n1.x - n2.x)
@inline normal_to_2Dline(I::NTuple{2,Integer},nodes::Vector{<:Vec2D}) = normal_to_2Dline(nodes[I[1]], nodes[I[2]])

@inline line_center(n1::AbstractVec,n2::AbstractVec) = 0.5*(n1+n2)
@inline line_center(I::NTuple{2,Integer},nodes::Vector{<:Vec2D}) = line_center(nodes[I[1]], nodes[I[2]])