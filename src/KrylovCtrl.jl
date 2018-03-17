
mutable struct KrylovCtrl{T,verb,min_iter,max_iter,abstol,reltol,bad_num,bad_val,kill_NaN}
  err0::T
  cnt::Int
  bad::Int
  exitval::String
end

function get_tiny(::Type{KrylovCtrl{T,verb,min_iter,max_iter,abstol,reltol,bad_num,bad_val,kill_NaN}}) where {T,verb,min_iter,max_iter,abstol,reltol,bad_num,bad_val,kill_NaN}
  if T === Float32
    tiny = Float32(1e-32)
  elseif T === Float64
    tiny = 1e-60
  end
  return tiny
end

@inline get_tiny(c::KrylovCtrl) = get_tiny(typeof(c))

for s in (:verb,:min_iter,:max_iter,:abstol,:reltol,:bad_num,:bad_val,:kill_NaN)
  f = Symbol(:get_,s)
  @eval begin
    ($f)(::Type{KrylovCtrl{T,verb,min_iter,max_iter,abstol,reltol,bad_num,bad_val,kill_NaN}}) where {T,verb,min_iter,max_iter,abstol,reltol,bad_num,bad_val,kill_NaN} = ($s)

    @inline ($f)(c::KrylovCtrl) = ($f)(typeof(c))
  end
end

function file_to_dict(i::AbstractString)
  open(i) do f
    d = Dict{Symbol,Any}()
    for l in eachline(f)
      a = split(l)
      if length(a) == 2
        d[Symbol(replace(a[1],":",""))] = a[2]
      elseif length(a) > 2
        d[Symbol(replace(a[1],":",""))] = (a[2:end]...)
      end
    end
    return d
  end
end

function KrylovCtrl{T}(infilename::AbstractString) where {T<:Real}
  d = file_to_dict(infilename)
  verb = parse(Int,d[:verbosity])
  min_iter = Int(parse(d[:min_max_iter][1]))
  max_iter = Int(parse(d[:min_max_iter][2]))
  abstol = parse(T,d[:abstol])
  reltol = parse(T,d[:reltol])
  badnum = parse(Int,d[:badtol][1])
  badval = parse(T,d[:badtol][2])
  kill_NaN = parse(Bool,d[:kill_NaN])

  return KrylovCtrl{T,verb,min_iter,max_iter,abstol,reltol,badnum,badval,kill_NaN}(0.0,0,0,"null")
end

KrylovCtrl(infile::AbstractString) = KrylovCtrl{Float64}(infile)
KrylovCtrl() = KrylovCtrl("PCG.set")

function init(ctrl::KrylovCtrl,ierr)
  ctrl.err0 = abs(ierr)
  ctrl.cnt = 0
  ctrl.bad = 0
  ctrl.exitval = "null"
  return nothing
end

function test(ctrl::KrylovCtrl,rho,beta=0.)
  ctrl.cnt += 1
  cnt = ctrl.cnt

  cnt <= get_min_iter(ctrl) && return true

  if beta > get_bad_val(ctrl)
    ctrl.bad += 1
  else
    ctrl.bad = 0
  end

  ar = abs(rho)

  if cnt > get_max_iter(ctrl)
    ctrl.exitval = "max"
    return false
  elseif ar < get_abstol(ctrl)
    ctrl.exitval = "abs"
    return false
  elseif ar < (ctrl.err0 * get_reltol(ctrl))
    ctrl.exitval = "rel"
    return false
  elseif ar < float(cnt)*get_tiny(ctrl)
    ctrl.exitval = "!!! round"
    return false
  elseif ctrl.bad >= get_bad_num(ctrl)
    ctrl.exitval = "!!! bad"
    return false
  end

  return true
end