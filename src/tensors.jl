@inline Base.@propagate_inbounds Base.getindex(a::AbstractTen,I::Integer) = getfield(a,I)
@inline Base.@propagate_inbounds Base.getindex(a::AbstractTen,I::Integer,J::Integer) = getfield(a, I + 2*(J-1))

struct Ten2D{T<:Number} <: AbstractTen{T}
    xx::T
    yx::T
    xy::T
    yy::T 
end

Base.size(a::Ten2D) = 
    (2,2)

Base.length(a::Ten2D) = 
    4

Base.linearindices(a::Ten2D) = 
    Base.OneTo(4)

Base.indices(a::Ten2D) = 
    (Base.OneTo(2),Base.OneTo(2))

Base.zero(a::Type{Ten2D{T}}) where {T} = 
    Ten2D{T}(zero(T),zero(T),zero(T),zero(T))

Ten2D(xx,yx,xy,yy) = 
    Ten2D(promote(xx,yx,xy,yy)...)

Base.:+(a::Ten2D{T},b::Ten2D{T2}) where {T,T2} = 
    Ten2D{promote_type(T,T2)}(a.xx+b.xx, a.yx+b.yx, a.xy+b.xy, a.yy+b.yy)

Base.:-(a::Ten2D{T},b::Ten2D{T2}) where {T,T2} = 
    Ten2D{promote_type(T,T2)}(a.xx-b.xx, a.yx-b.yx, a.xy-b.xy, a.yy-b.yy)

Base.:*(a::T,b::Ten2D{T2}) where {T<:Number,T2} = 
    Ten2D{promote_type(T,T2)}(a*b.xx, a*b.yx, a*b.xy, a*b.yy)

Base.:*(b::Ten2D{T},a::T2) where {T<:Number,T2} = 
    Ten2D{promote_type(T,T2)}(a*b.xx, a*b.yx, a*b.xy, a*b.yy)

Base.:/(b::Ten2D{T},a::T2) where {T<:Number,T2} = 
    Ten2D{promote_type(T,T2)}(b.xx/a, b.yx/a, b.xy/a, b.yy/a)

Base.dot(a::Ten2D,b::Ten2D) = 
    muladd(a.xx, b.xx, muladd(a.yx, b.yx, muladd(a.xy, b.xy, a.yy*b.yy)))

Base.norm(a::Ten2D) = 
    sqrt(2dot(a,a))

Base.dot(a::Ten2D{T},b::Vec2D{T2}) where {T,T2} = 
    Vec2D{promote_type(T,T2)}(a.xx*b.x + a.xy*b.y, a.yx*b.x + a.yy * b.y)

Base.dot(b::Vec2D{T},a::Ten2D{T2}) where {T,T2} = 
    Vec2D{promote_type(T,T2)}(a.xx*b.x + a.yx*b.y, a.xy*b.x + a.yy * b.y)

Base.:*(a::Vec2D{T},b::Vec2D{T2}) where {T,T2} = 
    Ten2D{promote_type(T,T2)}(a.x*b.x, a.y*b.x, a.x*b.y, a.y*b.y)