
area(face::AbstractFace2Node1D,nodes::Vector{<:Vec1D}) = 1 
area(face::AbstractFace2Node2D,nodes::Vector{<:Vec2D}) =  normal_to_2Dline(face.ind,nodes)

center(face::AbstractFace2Node1D,nodes::Vector{<:Vec1D}) = nodes[face[1]] 
center(face::AbstractFace2Node2D,nodes::Vector{<:Vec2D}) =  line_center(face.ind,nodes)

struct Face2Node{N} <: AbstractFace2Node{N}
  ind::NTuple{N,UInt}
end

abstract type AbstractFace2Cell <: NTupleWrap{2,UInt} end

abstract type AbstractBFace2Cell <: NTupleWrap{1,UInt} end

struct Face2Cell <: AbstractFace2Cell
  ind::Tuple{UInt,UInt}
end


struct BFace2Cell <: AbstractBFace2Cell
  ind::Tuple{UInt}
end

struct Face2CellLoop{NBF,NF,VecType}
  bf2c::Vector{BFace2Cell}
  bft::Vector{UInt} # List of of boundary faces types
  bfn::Vector{VecType} # List of normals of boundary faces
  bfc::Vector{VecType} # List of centers of boundary faces
  bcv::Vector{Float64} # List of volumes of boundary faces owner cells
  bccenter::Vector{VecType} # List of centers of boundary faces owner cells
  f2c::Vector{Face2Cell}
  fn::Vector{VecType} # List of normals of faces
  fc::Vector{VecType} # List of centers of faces
  cv::Vector{NTuple{2,Float64}} # List of volumes of faces owner cells
  ccenter::Vector{NTuple{2,VecType}} # List of centers of faces owner cells
end
