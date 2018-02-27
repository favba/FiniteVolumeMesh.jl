
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
