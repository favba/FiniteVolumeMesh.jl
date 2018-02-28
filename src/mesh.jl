
struct HomogeneousMesh{VecType<:AbstractVec,CellType<:AbstractCell2Node,CNF,FaceType<:AbstractFace2Node,NBF,NF} <: AbstractMesh
  nodes::Vector{VecType}
  cells::Vector{CellType} # c2n
  faces::Vector{FaceType} # f2n
  bfaces::Vector{FaceType} # bf2n
  #c2c::Vector{UInt}
  c2f::Vector{NTuple{CNF,Int}}
#  f2c::Vector{Face2Cell}
#  bf2c::Vector{BFace2Cell}
  f2cloops::Face2CellLoop{NBF,NF,VecType}
end

function HomogeneousMesh(inputfile::AbstractString)
  nodes, cells, bc_cell, bfaces, bc_face = readneutral(inputfile)
  faces, c2f, bf2c, bfn, bfc, bcv, bccenter = cell_connectivity(cells,bfaces,bc_face,nodes)
  f2c, fn, fc, cv, ccenter = face_connectivity(faces,c2f,cells,nodes)
  CellType = eltype(cells)
  CNF = nfaces(CellType)
  f2cloops = Face2CellLoop{length(bfaces),length(faces),eltype(nodes)}(bf2c,bc_face,bfn,bfc,bcv,bccenter,f2c,fn,fc,cv,ccenter)
  return HomogeneousMesh{eltype(nodes),CellType,CNF,eltype(faces),length(bfaces),length(faces)}(nodes,cells,faces,bfaces,c2f,f2cloops)
end

HomogeneousMesh(d::Dict) = HomogeneousMesh(d[:meshfilename])

function cell_connectivity(cells::Vector{<:TriangleCell},bfaces,bc_face,nodes)
  NC = length(cells)
  NBF = length(bfaces)

  faces = Vector{Face2Node{2}}(0)
  c2f = Vector{Tuple{Int,Int,Int}}(NC)
  bf2c = Vector{BFace2Cell}(NBF)
  bfn = Vector{Vec2D{Float64}}(NBF)
  bfc = Vector{Vec2D{Float64}}(NBF)
  bccenter = Vector{Vec2D{Float64}}(NBF)
  bcv = Vector{Float64}(NBF)
  
  nnbc = 0
  for (i,cell) in enumerate(cells)
    t1,t2,t3 = cell_faces(cell)

    t1notinfaces=true #future work: make this a tnotinfaces = zeros(Bool,NN)
    t2notinfaces=true
    t3notinfaces=true

    c2f_i1 = 0
    c2f_i2 = 0
    c2f_i3 = 0

    nnbc < NBF && for (j,el) in enumerate(bfaces)
      if t1 == el #future work make this a for t in tn, t == el
        c2f_i1 = -j
        bf2c[j] = BFace2Cell((i,))
        bfn[j] = area(t1,nodes) #future work: make this a call to area(face,nodes)
        bfc[j] = center(t1,nodes)
        bcv[j] = 1/volume(cell,nodes)
        bccenter[j] = center(cell,nodes)
        t1notinfaces = false
        nnbc+=1
      elseif t2 == el
        c2f_i2 = -j
        bf2c[j] = BFace2Cell((i,))
        bfn[j] = area(t2,nodes)
        bfc[j] = center(t2,nodes)
        bcv[j] = 1/volume(cell,nodes)
        bccenter[j] = center(cell,nodes)
        t2notinfaces = false
        nnbc+=1
      elseif t3 == el
        c2f_i3 = -j
        bf2c[j] = BFace2Cell((i,))
        bfn[j] =  area(t3,nodes)
        bfc[j] = center(t3,nodes)
        bcv[j] = 1/volume(cell,nodes)
        bccenter[j] = center(cell,nodes)
        t3notinfaces = false
        nnbc+=1
      end
    end
    
    for (j,el) in enumerate(faces)
      if (t1[2], t1[1]) == el
        c2f_i1 = j
        t1notinfaces = false
      elseif (t2[2], t2[1]) == el
        c2f_i2 = j
        t2notinfaces = false
      elseif (t3[2], t3[1]) == el
        c2f_i3 = j
        t3notinfaces = false
      end
    end

    if t1notinfaces
      push!(faces,t1)
      c2f_i1 = length(faces)
    end

    if t2notinfaces 
      push!(faces,t2)
      c2f_i2 = length(faces)
    end

    if t3notinfaces 
      push!(faces,t3)
      c2f_i3 = length(faces)
    end

    @assert (c2f_i1 != 0) & (c2f_i2 != 0) & (c2f_i3 != 0)
    c2f[i] = (c2f_i1, c2f_i2, c2f_i3)
  end
  return faces, c2f, bf2c, bfn, bfc, bcv, bccenter
end

function face_connectivity(faces::Vector{<:Face2Node},c2f,cells::Vector{T},nodes) where {T<:AbstractCell2Node}
  NF = length(faces)
  NC = length(c2f)
  #VecType = nodetype(T)
  f2c = Vector{Face2Cell}(NF)
  fn = Vector{Vec2D{Float64}}(NF)
  fc = Vector{Vec2D{Float64}}(NF)
  cv = Vector{NTuple{2,Float64}}(NF)
  ccenter = Vector{NTuple{2,Vec2D{Float64}}}(NF)
 
  for i=1:NF
    notfound = true
    j=0
    cell1 = UInt(0) 
    cell2 = UInt(0) 

    while notfound
      j+=1
      if i in c2f[j]
        notfound = false
        cell1 = UInt(j)
      end
    end
    
    for l=(j+1):NC
      if i in c2f[l]
        cell2 = UInt(l)
        break
      end
    end

    @assert (cell1 != 0) & (cell2 != 0)
    f2c[i] = Face2Cell((cell1, cell2))
    fn[i] =  area(faces[i],nodes) 
    fc[i] =  center(faces[i],nodes) 
    cv[i] =  (1/volume(cells[cell1],nodes), 1/volume(cells[cell2],nodes))
    ccenter[i] =  (center(cells[cell1],nodes), center(cells[cell2],nodes))

  end

  return f2c, fn, fc, cv, ccenter
end