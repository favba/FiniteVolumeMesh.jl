using FiniteVolumeMesh
using Base.Test

mesh = HomogeneousMesh("mesh.neu")

Tf = ones(length(mesh.faces))
Tbf = ones(length(mesh.bfaces))
∇T = FiniteVolumeMesh.gradient(Tf,Tbf,mesh)
@testset "Pure Gauss gradient 1st order accuracy" begin
  @test maximum(maximum,∇T) <= 1e-14
end

@testset "Pure Gauss gradient 2nd order accuracy" begin
  Tf .= (x->((0.5*(mesh.nodes[x[1]] + mesh.nodes[x[2]])).x)).(mesh.faces)
  Tbf .= (x->((0.5*(mesh.nodes[x[1]] + mesh.nodes[x[2]])).x)).(mesh.bfaces)
  FiniteVolumeMesh.gradient!(∇T,Tf,Tbf,mesh.f2cloops)
  @test minimum(x->x≈Vec2D(1.,0.),∇T)

  Tf .= (x->((0.5*(mesh.nodes[x[1]] + mesh.nodes[x[2]])).y)).(mesh.faces)
  Tbf .= (x->((0.5*(mesh.nodes[x[1]] + mesh.nodes[x[2]])).y)).(mesh.bfaces)
  FiniteVolumeMesh.gradient!(∇T,Tf,Tbf,mesh.f2cloops)
  @test minimum(x->x≈Vec2D(0.,1.),∇T)

  Tf .= (x->((0.5*(mesh.nodes[x[1]] + mesh.nodes[x[2]])).x)).(mesh.faces) .+ (x->((0.5*(mesh.nodes[x[1]] + mesh.nodes[x[2]])).y)).(mesh.faces) 
  Tbf .= (x->((0.5*(mesh.nodes[x[1]] + mesh.nodes[x[2]])).x)).(mesh.bfaces) .+ (x->((0.5*(mesh.nodes[x[1]] + mesh.nodes[x[2]])).y)).(mesh.bfaces) 
  FiniteVolumeMesh.gradient!(∇T,Tf,Tbf,mesh.f2cloops)
  @test minimum(x->x≈Vec2D(1.,1.),∇T)

end