
struct Homogeneous2DMesh{CellType} <: Abstract2DMesh
  nodes::Vector{Node2D{Float64}}
  cells::Vector{CellType} # c2n
  faces::Vector{Tuple{UInt,UInt}} # f2n
  bfaces::Vector{Tuple{UInt,UInt}} # bf2n
  #c2c::Vector{UInt}
  c2f::Vector{Tuple{Int,Int,Int}}
  f2c::Vector{Face{Node2D{Float64}}}
  bf2c::Vector{BoundaryFace{Node2D{Float64}}}
end

function Homogeneous2DMesh(inputfile::String)
  nodes, cells, bc_cell, bfaces, bc_face = readneutral(inputfile)
  faces, c2f, bf2c = cell_connectivity(cells,bfaces,bc_face,nodes)
  f2c = face_connectivity(faces,c2f,cells,nodes)
  return Homogeneous2DMesh{eltype(cells)}(nodes,cells,faces,bfaces,c2f,f2c,bf2c)
end

function cell_connectivity(cells::Vector{TriangleCell},bfaces,bc_face,nodes)
  NC = length(cells)
  NBF = length(bfaces)

  faces = Vector{Tuple{UInt,UInt}}(0)
  c2f = Vector{Tuple{Int,Int,Int}}(NC)
  bf2c = Vector{BoundaryFace{Node2D{Float64}}}(NBF)
  
  nnbc = 0
  for (i,cell) in enumerate(cells)
    t1,t2,t3 = cell_faces(cell)

    t1notinfaces=true
    t2notinfaces=true
    t3notinfaces=true

    c2f_i1 = 0
    c2f_i2 = 0
    c2f_i3 = 0

    nnbc < NBF && for (j,el) in enumerate(bfaces)
      if t1 == el
        c2f_i1 = -j
        bf2c[j] = BoundaryFace{Node2D{Float64}}(i, normal_to_2Dline(t1,nodes), volume(cell,nodes),bc_face[j])
        t1notinfaces = false
        nnbc+=1
      elseif t2 == el
        c2f_i2 = -j
        bf2c[j] = BoundaryFace{Node2D{Float64}}(i, normal_to_2Dline(t2,nodes), volume(cell,nodes),bc_face[j])
        t2notinfaces = false
        nnbc+=1
      elseif t3 == el
        c2f_i3 = -j
        bf2c[j] = BoundaryFace{Node2D{Float64}}(i, normal_to_2Dline(t3,nodes), volume(cell,nodes),bc_face[j])
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
  return faces, c2f, bf2c
end

function face_connectivity(faces::Vector{Tuple{UInt,UInt}},c2f,cells::Vector{TriangleCell},nodes)
  NF = length(faces)
  NC = length(c2f)
  f2c = Vector{Face{Node2D{Float64}}}(NF)

  for i=1:NF
    notfound = true
    j=0
    cell1 = UInt(0) 
    cell2 = UInt(0) 
    while notfound
      j+=1
      ce = c2f[j]

      if ce[1] == i
        notfound = false
        cell1 = UInt(j)
      elseif ce[2] == i
        notfound = false
        cell1 = UInt(j)
      elseif ce[3] == i
        notfound = false
        cell1 = UInt(j)
      end

    end
    
    for l=(j+1):NC
      ce = c2f[l]

      if ce[1] == i
        cell2 = UInt(l)
        break
      elseif ce[2] == i
        cell2 = UInt(l)
        break
      elseif ce[3] == i
        cell2 = UInt(l)
        break
      end

    end

    @assert (cell1 != 0) & (cell2 != 0)
    f2c[i] = Face{Node2D{Float64}}(cell1, cell2,  
      normal_to_2Dline(faces[i],nodes), 
      volume(cells[cell1],nodes), volume(cells[cell2],nodes))

  end

  return f2c
end