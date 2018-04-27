abstract type AbstractProblem end

struct CellProblemAdvecTemp{MeshType,TbfType,TypeK,TypeS,VType,BcondType,LaplacianStruct,AdvectionStruct} <: AbstractProblem
    Tc::Vector{Float64}
    k::TypeK
    s::TypeS
    ρC::Float64
    u::VType
    bcond::BcondType
    mesh::MeshType
    laplacian!::LaplacianStruct
    advection!::AdvectionStruct
    ∇Tc::Vector{Vec2D{Float64}}
    Tbf::TbfType
end

function CellProblemAdvecTemp(Tc,mesh,d)
    bcond = boundary_conditions(d)

    k = ConstVec{d[:conductivity]}()
    s = ConstVec{d[:source]}()
    ρC = d[Symbol("rho*C")]
    uvec = ConstVec{d[:velocity]}()
    ∇Tc = zeros(Vec2D{Float64},length(Tc))
    Tbf = FieldAtBoundary(Tc,mesh,bcond)

    advection = get_advection_method(d)
    laplacian = get_laplacian_method(Tc,bcond,k,mesh,d)
    types = typeof.((mesh,Tbf,k,s,uvec,bcond,laplacian,advection))
    return CellProblemAdvecTemp{types...}(Tc,k,s,ρC,uvec,bcond,mesh,laplacian,advection,∇Tc,Tbf)
end

@inline needs_gradient_calculation(p) = (typeof(p.advection!) === UpWind2ndOrder) && (typeof(p.laplacian!) === MethodB) && (typeof(p.u) <: ConstVec && p.u[1] == zero(eltype(p.u)))

function calculate_rhs!(rhs,p::CellProblemAdvecTemp)
    fill!(rhs,zero(eltype(rhs)))

    is_implicit(p.laplacian!) || p.laplacian!(rhs,p)

    needs_gradient_calculation(p) && gradient!(p.∇Tc,FaceSimpleInterpolation(p.Tc,p.mesh),p.Tbf,p.mesh.f2cloops)

    if !(typeof(p.u) <: ConstVec && p.u[1] == zero(eltype(p.u))) 
        p.advection!(rhs,p)
    end
    if !(typeof(p.s) <: ConstVec && p.s[1] == zero(eltype(p.s))) 
        add_source!(rhs,p)
    end

    is_implicit(p.laplacian!) && p.laplacian!(rhs,p)

    return nothing
end

function add_source!(rhs,p)
    s = p.s
    @inbounds @simd for i=1:length(rhs)
        rhs[i] += s[i]
    end
end

struct StokesProblem{UfType,UbfType,MeshType,uBcondType,LaplacianStruct,Ma,PMa,PMap} <: AbstractProblem
    u::Vec2DArray{Float64}
    u_old::Vec2DArray{Float64}
    uf::UfType
    ubf::UbfType
    δu::Vec2DArray{Float64}
    ru::Vec2DArray{Float64}
    p::Vector{Float64}
    δp::Vector{Float64}
    ν::ConstVec{Float64}
    ρ::ConstVec{Float64}
    bcond::uBcondType
    mesh::MeshType
    laplacian!::LaplacianStruct
    ∇u::Vector{FiniteVolumeMesh.Ten2D{Float64}}
    A::Ma
    pcgA::PMa
    s::Vector{Float64}
    Ap::PoissonP{MeshType}
    pcgAp::PMap
    dt::Float64
end

function StokesProblem(u,mesh,d)
    u_old = similar(u)
    uf = FaceSimpleInterpolation(u,mesh)
    bcond = uboundary_conditions(d)
    ubf = FieldAtBoundary(u,mesh,bcond)
    δu = similar(u)
    ru = similar(u)
    ν = ConstVec(d[:viscosity])
    ρ = ConstVec(d[:density])
    ∇u = zeros(Ten2D{Float64},length(u))
    p = zeros(length(u))
    δp = zeros(length(u))

    laplacian = get_laplacian_method(u,bcond,ν,mesh,d)
    A = aIpDbG(d,mesh)
    pcgA = PCG(δu)
    pcgAp = PCG(δp)
    types = typeof.((uf,ubf,mesh,bcond,laplacian,A,pcgA,pcgAp))

    s = zeros(length(mesh.cells))
    Ap = PoissonP(d,mesh)
    dt = d[:dt]
    return StokesProblem{types...}(u,u_old,uf,ubf,δu,ru,p,δp,ν,ρ,bcond,mesh,laplacian,∇u,A,pcgA,s,Ap,pcgAp,dt)
end

struct NSProblem{UfType,UbfType,MeshType,uBcondType,LaplacianStruct,AdvecStruct,Ma,PMa,PMap,TBcondType,TLaplacianStruct,TAdvectionStruct,TbfType,MaT,PMaT} <: AbstractProblem
    u::Vec2DArray{Float64}
    uf::UfType
    ubf::UbfType
    δu::Vec2DArray{Float64}
    ru::Vec2DArray{Float64}
    p::Vector{Float64}
    δp::Vector{Float64}
    ν::ConstVec{Float64}
    ρ::ConstVec{Float64}
    bcond::uBcondType
    mesh::MeshType
    laplacian!::LaplacianStruct
    uadvection!::AdvecStruct
    ∇u::Vector{FiniteVolumeMesh.Ten2D{Float64}}
    A::Ma
    pcgA::PMa
    sp::Vector{Float64}
    Ap::PoissonP{MeshType}
    pcgAp::PMap
    dt::Float64
    Tc::Vector{Float64}
    ∇Tc::Vec2DArray{Float64}
    rT::Vector{Float64}
    δT::Vector{Float64}
    k::ConstVec{Float64}
    s::ConstVec{Float64}
    ρC::Float64
    Tbcond::TBcondType
    Tlaplacian!::TLaplacianStruct
    Tadvection!::TAdvectionStruct
    Tbf::TbfType
    AT::MaT
    pcgAT::PMaT
end

function NSProblem(u,mesh,d)
    uf = FaceSimpleInterpolation(u,mesh)
    bcond = uboundary_conditions(d)
    ubf = FieldAtBoundary(u,mesh,bcond)
    δu = Vec2DArray{Float64}(length(u))
    ru = Vec2DArray{Float64}(length(u))
    ν = ConstVec(d[:viscosity])
    ρ = ConstVec(d[:density])
    ∇u = zeros(Ten2D{Float64},length(u))
    p = zeros(length(u))
    δp = zeros(length(u))

    laplacian = get_laplacian_method(u,bcond,ConstVec{Float64}(0.5*ν[1]),mesh,d)
    A = aIpDbG(ConstVec(d[:density]/d[:dt]),-0.5*ν[1],d,mesh)
    pcgA = PCG(δu,KrylovCtrl("vPCG.set"))
    pcgAp = PCG(δp,KrylovCtrl("pPCG.set"))
    uadvection = UpWindAdvection(u,ubf,u,ubf,mesh)

    sp = zeros(length(mesh.cells))
    Ap = PoissonP(d,mesh)
    dt = d[:dt]

    Tc = zeros(Float64,length(mesh.cells))
    ∇Tc = Vec2DArray{Float64}(length(u))
    rT = zeros(Float64,length(mesh.cells))
    dT = zeros(Float64,length(mesh.cells))
    Tbcond = boundary_conditions(d)

    k = ConstVec{Float64}(d[:conductivity])
    s = ConstVec{Float64}(d[:source])
    ρC = d[Symbol("rho*C")]
    Tbf = FieldAtBoundary(Tc,mesh,Tbcond)

    Tadvection = UpWindAdvection(Tc,Tbf,u,ubf,mesh)
    Tlaplacian = get_laplacian_method(Tc,Tbcond,ConstVec{Float64}(0.5*k[1]),mesh,d)
    AT = aIpDbG(ConstVec{Float64}(ρC/dt),-0.5*k[1],d,mesh)
    pcgAT = PCG(dT,KrylovCtrl("tPCG.set"))

    types = typeof.((uf,ubf,mesh,bcond,laplacian,uadvection,A,pcgA,pcgAp,Tbcond,Tlaplacian,Tadvection,Tbf,AT,pcgAT))
    return NSProblem{types...}(u,uf,ubf,δu,ru,p,δp,ν,ρ,bcond,mesh,laplacian,uadvection,∇u,A,pcgA,sp,Ap,pcgAp,dt,Tc,∇Tc,rT,dT,k,s,ρC,Tbcond,Tlaplacian,Tadvection,Tbf,AT,pcgAT)
end