__precompile__()
module FiniteVolumeMesh

export HomogeneousMesh, Vec, Vec2D, TriangleCell, read_input, CellProblemAdvecTemp, PCG

include("utils.jl")

abstract type AbstractVec{T} <: DenseArray{T,1} end
abstract type AbstractTen{T} <: DenseArray{T,2} end

abstract type AbstractCell2Node{NN,NF} <: NTupleWrap{NN,UInt} end
abstract type AbstractCell2Node2D{NN,NF} <: AbstractCell2Node{NN,NF} end

abstract type AbstractFace2Node{N} <: NTupleWrap{N,UInt} end
const AbstractFace2Node1D = AbstractFace2Node{1}
const AbstractFace2Node2D = AbstractFace2Node{2}

abstract type AbstractMesh end

include("KrylovCtrl.jl")
include("pcg.jl")
include("vec.jl")
include("vecarray.jl")
include("tensors.jl")
include("tensorsarray.jl")
include("fieldtypes.jl")
include("geometry.jl")
include("cells.jl")
include("faces.jl")
include("boundaryconditions.jl")
include("mesh.jl")
include("MatrixLike.jl")
include("advection.jl")
include("input.jl")
include("readneutral.jl")
include("writevtk.jl")
include("derivatives.jl")
include("solvers.jl")
include("laplacian.jl")

end # module
