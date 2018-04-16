using Base.RefValue
struct PCG{A,CtrlType}
    ctrl::CtrlType
    pp::A
    zz::A
    ww::A
    δ::RefValue{Float64}
    α::RefValue{Float64}
    β::RefValue{Float64}
    η::RefValue{Float64}
    ηm::RefValue{Float64}
    err::RefValue{Float64}
end

function PCG(xx,ictrl)
    pp = similar(xx)
    zz = similar(xx)
    ww = similar(xx)
    return PCG{typeof(xx),typeof(ictrl)}(ictrl,pp,zz,ww,Ref(0.),Ref(0.),Ref(0.),Ref(0.),Ref(0.),Ref(0.))
end

function PCG(xx)
    ictrl = KrylovCtrl()
    return PCG(xx,ictrl)
end

function solve!(pcg::PCG,Mat,xx,rr,ini::Bool=false)
    header(pcg,Mat,xx,rr,ini)
    out = 0
    while (test(pcg.ctrl,pcg.err[],pcg.β[]) && out == 0)
        out = iter(pcg,Mat,xx,rr)
        output(pcg)
    end
    out != 0 && println("Exite jumped out at ",out)
    footer(pcg)
end

function header(pcg,Mat,xx,rr,ini::Bool)
    if ini
        pcg.η[] = evalT(Mat,pcg.zz,xx)
        @inbounds for i in linearindices(xx)
            rr[i] -= pcg.zz[i]
        end
    else
        @inbounds for i in linearindices(xx)
            xx[i] = zero(eltype(xx))
        end
    end
    pcg.η[] = evalP(Mat,pcg.zz,rr)
    pp = pcg.pp
    zz = pcg.zz
    @inbounds for i in linearindices(xx)
        pp[i] = zz[i]
    end
    pcg.err[] = pcg.η[]
    init(pcg.ctrl,pcg.err[])
    output(pcg)
end

function iter(pcg::PCG,Mat,xx,rr)
    pcg.δ[] = evalT(Mat,pcg.ww,pcg.pp)
    δ = pcg.δ[]
    tiny = get_tiny(pcg.ctrl)
  
    abs(δ) < tiny && return 1
  
    pcg.α[] = pcg.η[]/δ
    α = pcg.α[]
    pp = pcg.pp
    ww = pcg.ww
    @inbounds for i in linearindices(xx)
        xx[i] += α*pp[i]
        rr[i] -= α*ww[i]
    end
    pcg.ηm[] = pcg.η[]
    ηm = pcg.ηm[]
    pcg.η[]  = evalP(Mat,pcg.zz,rr)
    η = pcg.η[]
    abs(η) < tiny && return 2

    pcg.β[] = η/ηm
    β = pcg.β[]
    zz = pcg.zz
    for i in linearindices(xx)
        pp[i] = zz[i] + β*pp[i]
    end

    pcg.err[] = η
    return 0
end

function footer(pcg)
    verb = get_verb(pcg.ctrl)
    if verb >= 1
        println("PCG iterations = ",pcg.ctrl.cnt - 1)
        println(pcg.ctrl.exitval, "\t exit. Final err = ",pcg.err[]," ee = ", sqrt(abs(pcg.err[]/pcg.ctrl.err0)))
        println()
    end
end

function output(pcg::PCG)
    ctrl = pcg.ctrl
    verb = get_verb(ctrl)
    cnt = ctrl.cnt
    if cnt == 0
        if verb >= 1
            println()
            println("Initial err =",ctrl.err0)
            println("atol = ",get_abstol(ctrl),"\t","rtol = ",get_reltol(ctrl))
            verb >= 2 && println("  m,    err,    alpha,    beta,    delta  ")
        end
    else
        verb >= 2 && println(join((cnt,pcg.err[],pcg.α[],pcg.β[],pcg.δ[]),"\t")) 
    end
end

function evalT(A::AbstractArray{T,2},out::AbstractVector{T2},input::AbstractVector{T2}) where {T,T2}
    A_mul_B!(out,A,input)
    δ = out⋅input
    return δ
end

function evalP(A::AbstractArray{T,2},out::AbstractVector{T2},input::AbstractVector{T2}) where {T,T2}
    δ = zero(T2)
    for i in linearindices(input)
        out[i] = input[i]/A[i,i]
        δ += out[i]*input[i]
    end
    return δ
end
