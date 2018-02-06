function writevtk(outputname::String,nodes::Vector{<:AbstractNode},cells::Vector{<:TriangleCell})
  info("Writing to $outputname")
  open(outputname,"w") do output
    NN = length(nodes)
    NC = length(cells)
    print(output,"# vtk DataFile Version 3.0\nThis file was generated from a neutral file\nASCII\n")
    print(output,"DATASET UNSTRUCTURED_GRID\n")
    print(output,"POINTS ",NN," float\n")
    for node in nodes
      println(output,xpos(node)," ",ypos(node)," ",zpos(node))
    end
    println(output,"CELLS ",NC," ",4*NC)
    for cell in cells
      #Juila starts indexing at 1, so we subtract 1 from node index.
      println(output,"3 ",cell.node1-1, " ", cell.node2-1, " ",cell.node3-1)
    end
    println(output,"CELL_TYPES ",NC)
    for cell in cells
      println(output,vtkcelltype(cell))
    end
  end
end

function append_points_scalar_data(vtkfile::String,dataname::String,data::AbstractArray,header::Bool)
  open(vtkfile,"a") do out
    header && println(out,"POINT_DATA ",length(data))
    println(out,"SCALARS ",dataname," float 1")
    println(out,"LOOKUP_TABLE default")
    for el in data
      println(out,el)
    end
  end
end

function append_cells_scalar_data(vtkfile::String,dataname::String,data::AbstractArray,header::Bool)
  open(vtkfile,"a") do out
    header && println(out,"CELL_DATA ",length(data))
    println(out,"SCALARS ",dataname," float 1")
    println(out,"LOOKUP_TABLE default")
    for el in data
      println(out,el)
    end
  end
end