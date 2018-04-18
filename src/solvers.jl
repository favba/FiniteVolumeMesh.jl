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
    u::Vector{Vec2D{Float64}}
    uf::UfType
    ubf::UbfType
    δu::Vector{Vec2D{Float64}}
    ru::Vector{Vec2D{Float64}}
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
    uf = FaceSimpleInterpolation(u,mesh)
    bcond = uboundary_conditions(d)
    ubf = FieldAtBoundary(u,mesh,bcond)
    δu = similar(u)
    ru = similar(u)
    ν = ConstVec(d[:conductivity])
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
    return StokesProblem{types...}(u,uf,ubf,δu,ru,p,δp,ν,ρ,bcond,mesh,laplacian,∇u,A,pcgA,s,Ap,pcgAp,dt)
end

struct NSProblem{UfType,UbfType,MeshType,uBcondType,LaplacianStruct,AdvecStruct,Ma,PMa,PMap} <: AbstractProblem
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
    s::Vector{Float64}
    Ap::PoissonP{MeshType}
    pcgAp::PMap
    dt::Float64
end

function NSProblem(u,mesh,d)
    uf = FaceSimpleInterpolation(u,mesh)
    bcond = uboundary_conditions(d)
    ubf = FieldAtBoundary(u,mesh,bcond)
    δu = Vec2DArray{Float64}(length(u))
    ru = Vec2DArray{Float64}(length(u))
    ν = ConstVec(d[:conductivity])
    ρ = ConstVec(d[:density])
    ∇u = zeros(Ten2D{Float64},length(u))
    p = zeros(length(u))
    δp = zeros(length(u))

    laplacian = get_laplacian_method(u,bcond,ν,mesh,d)
    A = aIpDbG(d,mesh)
    pcgA = PCG(δu,KrylovCtrl("vPCG.set"))
    pcgAp = PCG(δp,KrylovCtrl("pPCG.set"))
    uadvection = UpWindAdvection(u,ubf,uf,ubf,mesh)
    types = typeof.((uf,ubf,mesh,bcond,laplacian,uadvection,A,pcgA,pcgAp))

    s = zeros(length(mesh.cells))
    Ap = PoissonP(d,mesh)
    dt = d[:dt]
    return NSProblem{types...}(u,uf,ubf,δu,ru,p,δp,ν,ρ,bcond,mesh,laplacian,uadvection,∇u,A,pcgA,s,Ap,pcgAp,dt)
end