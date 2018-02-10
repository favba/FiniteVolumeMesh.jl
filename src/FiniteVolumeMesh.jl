__precompile__()
module FiniteVolumeMesh

abstract type AbstractVec{T} <: DenseArray{T,1} end

abstract type AbstractCell end
abstract type Abstract2DCell <: AbstractCell end

abstract type AbstractFace end
abstract type AbstractBoundaryFace end

abstract type AbstractMesh end
abstract type Abstract2DMesh <: AbstractMesh end

include("vec.jl")
include("geometry.jl")
include("cells.jl")
include("faces.jl")
include("mesh.jl")
include("readneutral.jl")
include("writevtk.jl")

end # module
