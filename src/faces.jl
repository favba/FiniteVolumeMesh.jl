#=
Idea: 

AbstractFace2Node
AbstractBFace2Node
AbstractFace2Cell
AbstractBFace2Cell

struct Face2Node2D
struct Face2Cell{NodeType} (Face{NodeType} below)
struct BFace2Cell{NodeType} (BoundaryFace{NodeType} below)

struct Face2Node{PolygonType}
=#

struct Face32{NodeType} <: AbstractFace
  cell1::UInt32 # Owner cell index
  cell2::UInt32 # Neighbor cell index
  n::NodeType # Normal vector
  vc1::Float32 # Cell 1 Volume
  vc2::Float32 # Cell 2 Volume
end

struct Face{NodeType} <: AbstractFace
  cell1::UInt64 # Owner cell index
  cell2::UInt64 # Neighbor cell index
  n::NodeType # Normal vector
  vc1::Float64 # Cell 1 Volume
  vc2::Float64 # Cell 2 Volume
end

struct BoundaryFace{NodeType} <: AbstractBoundaryFace
  cell::UInt64 # Owner cell index
  n::NodeType # Normal vector
  vc::Float64 # Owner Cell Volume
  btype::UInt64 # Boundary type
end

area(face::Tuple{UInt,UInt},nodes::Vector{<:Vec2D}) =  normal_to_2Dline(face,nodes)