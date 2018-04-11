function read_input(filename::String)
  d =  Dict{Symbol,Any}()
  open(filename) do f
    # read meshfilename and outfilename options
    for i=1:2
      r = split(readline(f))
      d[Symbol(r[1])] = r[2]
    end
    # read maxsteps and printskip options
    for i=1:2
      r = split(readline(f))
      d[Symbol(r[1])] = parse(Int,r[2])
    end
    # read pause(1=true) option
    r = split(readline(f))
    d[Symbol(r[1])] = Bool(parse(Int,r[2]))
  
    #skip blank line 
    readline(f)

    #read dt and source
    for i=1:2
      r = split(readline(f))
      d[Symbol(r[1])] = parse(Float64,r[2])
    end

    #read gradtype
    r = split(readline(f)) 
    d[Symbol(r[1])] = r[2]

    #skip blank line 
    readline(f)

    #read condudctivity and rho*C
    for i=1:2
      r = split(readline(f))
      d[Symbol(r[1])] = parse(Float64,r[2])
    end
  
    #read velocity
    r = split(readline(f))
    d[Symbol(r[1])] = Vec2D(parse(Float64,r[2]), parse(Float64,r[3])) 

    #read viscosity and density
    for i=1:2
      r = split(readline(f))
      d[Symbol(r[1])] = parse(Float64,r[2])
    end
 
    #skip blank line
    readline(f)

    #read velocity
    r = split(readline(f))
    d[Symbol(r[1])] = parse(Bool,r[2])
    #read boundary conditions
    r = split(readline(f)) 
    r[1] == "BCs" || error("Error reading $filename file line 15. Expected parameter \"BCs\", got \"$(r[1])\" instead.")
    nbcs = parse(UInt,r[2])
    bcs = ()
    for i=1:nbcs
      r = split(readline(f))
      bc = Symbol(r[1])=>Vec2D(parse(Float64,r[2]), parse(Float64,r[3])) 
      bcs = (bcs...,bc)
    end
    d[:uBCs] = bcs
    #skip blank line
    readline(f)

    #read temperature
    r = split(readline(f))
    d[Symbol(r[1])] = parse(Bool,r[2])
    #read boundary conditions
    r = split(readline(f)) 
    r[1] == "BCs" || error("Error reading $filename file line 15. Expected parameter \"BCs\", got \"$(r[1])\" instead.")
    nbcs = parse(UInt,r[2])
    bcs2 = ()
    for i=1:nbcs
      r = split(readline(f))
      bc2 = Symbol(r[1])=>parse(Float64,r[2])
      bcs2 = (bcs2...,bc2)
    end
    d[:BCs] = bcs2
  end
  return d
end

read_input() = read_input("input.txt")

function show_parameters(d::Dict)
  msg = "Parameters from input:\n"
  for el in d
    msg *= string(el[1],": ",el[2],"\n")
  end
  msg*="\n"
  return msg
end

const available_laplacians = ("Method-B","Method-Gauss","Corr-Gauss","Mimetic-Diff","Implicit-Diff")

function get_laplacian_method(Tc,bcond,k,mesh,d)
  method = d[:gradtype]
  if method == "Method-B"
    return MethodB()
  elseif method == "Method-Gauss"
    return Gauss(Tc,mesh)
  elseif method == "Corr-Gauss"
    return CorrectedGauss(Tc,mesh)
  elseif method == "Mimetic-Diff" || method == "Mimetic_Diff"
    return CellMimetic(k,bcond,mesh)
  elseif method == "Implicit-Diff" || method == "Implicit_Diff"
    return ImplicitDiff(k,bcond,mesh,d)
  else
    warn("Gradtype \"",method,"\" not implemented. Using \"Corr-Gauss\". Available options are: ",available_laplacians)
    return CorrectedGauss(Tc,mesh)
  end
end

const available_advections = ("FullUpWind","UpWind2ndOrder")

function get_advection_method(d)
  if haskey(d,:advectype)
    method = d[:advectype]
    if method == "FullUpWind"
      return FullUpWind()
    elseif method == "UpWind2ndOrder"
      return UpWind2ndOrder()
    else
      warn("Advectype \"",method,"\" not implemented. Using \"UpWind2ndOrder\". Available options are: ",available_advections)
      return FullUpWind()
    end
  else
    return UpWind2ndOrder()
  end
end
