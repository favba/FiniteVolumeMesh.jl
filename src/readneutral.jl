
string_to_numbers(T::DataType,str::String) = parse.(T,split(str," ",keep=false))

function readneutral(inputfile::String)
  info("Reading $inputfile file")
  open(inputfile) do input

    NN = parse(Int,readline(input))
    nodes = Vector{Node2D{Float64}}(NN) #Pre-allocation a vector of Node2D{Float64} with NN elemetns
    
    for i=1:NN
      nodes[i] = Node2D(string_to_numbers(Float64,readline(input))...)
    end

    NC = parse(Int,readline(input))

    cells = Vector{TriangleCell}(NC)
    bc_cell = Vector{UInt}(NC)

    for i=1:NC
      a = string_to_numbers(UInt,readline(input)) 
      bc_cell[i] = a[1]
      cells[i] = TriangleCell(a[2],a[3],a[4])
    end

    NBF = parse(Int,readline(input))

    bf2n =  Vector{Tuple{UInt,UInt}}(NBF)
    bc_face = Vector{UInt}(NBF)

    for i=1:NBF
      a = string_to_numbers(UInt,readline(input))
      bc_face[i] = a[1]
      bf2n[i] = (a[2],a[3])
    end

    return nodes, cells, bc_cell, bf2n, bc_face
  end
end
