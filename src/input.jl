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

    #skip blank line
    readline(f)

    #read boundary conditions
    r = split(readline(f)) 
    r[1] == "BCs" || error("Error reading $filename file line 15. Expected parameter \"BCs\", got \"$(r[1])\" instead.")
    nbcs = parse(UInt,r[2])
    bcs = ()
    for i=1:nbcs
      r = split(readline(f))
      bc = Symbol(r[1])=>parse(Float64,r[2])
      bcs = (bcs...,bc)
    end
    d[:BCs] = bcs
  end
  return d
end

read_input() = read_input("input.txt")