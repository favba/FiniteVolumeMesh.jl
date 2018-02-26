@inline Base.@propagate_inbounds Base.getindex(a::AbstractVec,I::Integer) = getfield(a,I)
Base.IndexStyle(a::Type{<:AbstractVec}) = Base.IndexLinear()

@inline xpos(a::AbstractVec) = a.x
@inline ypos(a::AbstractVec) = a.y
@inline zpos(a::AbstractVec) = a.z

#pos of a node returns a tuple with the positions
pos(a::AbstractVec) = (xpos(a),ypos(a),zpos(a))

Base.norm(a::AbstractVec) = sqrt(muladd(xpos(a), xpos(a), muladd(ypos(a), ypos(a), zpos(a)^2)))
distance(a::AbstractVec,b::AbstractVec) = sqrt((xpos(b)-xpos(a))^2 + (ypos(b)-ypos(a))^2 + (zpos(b)-zpos(a))^2)

Base.:+(a::AbstractVec,b::AbstractVec) = Vec(xpos(a)+xpos(b), ypos(a)+ypos(b), zpos(a)+zpos(b))
Base.:-(a::AbstractVec,b::AbstractVec) = Vec(xpos(a)-xpos(b), ypos(a)-ypos(b), zpos(a)-zpos(b))
Base.:*(a::Number,b::AbstractVec) = Vec(a*xpos(b), a*ypos(b), a*zpos(b))
Base.:*(b::AbstractVec,a::Number) = Vec(a*xpos(b), a*ypos(b), a*zpos(b))
Base.:/(b::AbstractVec,a::Number) = Vec(xpos(b)/a, ypos(b)/a, zpos(b)/a)

Base.dot(a::AbstractVec,b::AbstractVec) = muladd(xpos(a), xpos(b), muladd(ypos(a), ypos(b), zpos(a)*zpos(b)))

Base.cross(a::AbstractVec,b::AbstractVec) = Vec(ypos(a)*zpos(b) - zpos(a)*ypos(b), zpos(a)*xpos(b) - xpos(a)*zpos(b), xpos(a)*ypos(b) - ypos(a)*xpos(b))

struct Vec{T<:Number} <: AbstractVec{T}
  x::T
  y::T
  z::T
end

Base.size(a::Vec) = (3,)
Base.length(a::Vec) = 3
Base.linearindices(a::Vec) = Base.OneTo(3)
Base.indices(a::Vec) = (Base.OneTo(3),)
Base.zero(a::Type{Vec{T}}) where {T} = Vec{T}(zero(T),zero(T),zero(T))
Vec(x,y,z) = Vec(promote(x,y,z)...)

struct Vec2D{T<:Number} <: AbstractVec{T}
  x::T
  y::T 
end

Base.size(a::Vec2D) = (2,)
Base.length(a::Vec2D) = 2
Base.linearindices(a::Vec2D) = Base.OneTo(2)
Base.indices(a::Vec2D) = (Base.OneTo(2),)
Base.zero(a::Type{Vec2D{T}}) where {T} = Vec2D{T}(zero(T),zero(T))

Vec2D(x,y) = Vec2D(promote(x,y)...)

@inline zpos(a::Vec2D{T}) where {T} = zero(T)

Base.:+(a::Vec2D{T},b::Vec2D{T2}) where {T,T2} = Vec2D{promote_type(T,T2)}(xpos(a)+xpos(b), ypos(a)+ypos(b))
Base.:-(a::Vec2D{T},b::Vec2D{T2}) where {T,T2} = Vec2D{promote_type(T,T2)}(xpos(a)-xpos(b), ypos(a)-ypos(b))
Base.:*(a::T,b::Vec2D{T2}) where {T<:Number,T2} = Vec2D{promote_type(T,T2)}(a*xpos(b), a*ypos(b))
Base.:*(b::Vec2D{T},a::T2) where {T,T2<:Number} = Vec2D{promote_type(T,T2)}(a*xpos(b), a*ypos(b))
Base.:/(b::Vec2D{T},a::T2) where {T,T2<:Number} = Vec2D{promote_type(T,T2)}(xpos(b)/a, ypos(b)/a)


Base.dot(a::Vec2D,b::Vec2D) = muladd(xpos(a), xpos(b), ypos(a)*ypos(b))
Base.norm(a::Vec2D) = sqrt(muladd(xpos(a), xpos(a), ypos(a)^2))
Base.cross(a::Vec2D,b::Vec2D) = xpos(a)*ypos(b) - ypos(a)*xpos(b)

struct Vec1D{T<:Number} <: AbstractVec{T}
  x::T
end

Base.size(a::Vec1D) = (1,)
Base.length(a::Vec1D) = 1
Base.linearindices(a::Vec1D) = Base.OneTo(1)
Base.indices(a::Vec1D) = (Base.OneTo(1),)
Base.zero(a::Type{Vec1D{T}}) where {T} = Vec1D{T}(zero(T))

@inline ypos(a::Vec1D{T}) where {T} = zero(T)
@inline zpos(a::Vec1D{T}) where {T} = zero(T)

Base.:+(a::Vec1D{T},b::Vec1D{T2}) where {T,T2} = Vec1D{promote_type(T,T2)}(xpos(a)+xpos(b))
Base.:-(a::Vec1D{T},b::Vec1D{T2}) where {T,T2} = Vec1D{promote_type(T,T2)}(xpos(a)-xpos(b))
Base.:*(a::T,b::Vec1D{T2}) where {T<:Number,T2} = Vec1D{promote_type(T,T2)}(a*xpos(b))
Base.:*(b::Vec1D{T},a::T2) where {T,T2<:Number} = Vec1D{promote_type(T,T2)}(a*xpos(b))
Base.:/(b::Vec1D{T},a::T2) where {T,T2<:Number} = Vec1D{promote_type(T,T2)}(xpos(b)/a)


Base.dot(a::Vec1D,b::Vec1D) = xpos(a)*xpos(b)
Base.norm(a::Vec1D) = abs(xpos(a))
Base.cross(a::Vec1D,b::Vec1D) = 0
