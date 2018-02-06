#This is a short way of defining one-line function, Ex: f(x) = x^2
#The function accepts any subtype(or child) of the AbstractNode type, which was defined above.
xpos(a::AbstractNode) = a.x
ypos(a::AbstractNode) = a.y
zpos(a::AbstractNode) = a.z

#pos of a node returns a tuple with the positions
pos(a::AbstractNode) = (xpos(a),ypos(a),zpos(a))

Base.norm(a::AbstractNode) = sqrt(muladd(xpos(a), xpos(a), muladd(ypos(a), ypos(a), zpos(a)^2)))
distance(a::AbstractNode,b::AbstractNode) = sqrt((xpos(b)-xpos(a))^2 + (ypos(b)-ypos(a))^2 + (zpos(b)-zpos(a))^2)

Base.:+(a::AbstractNode,b::AbstractNode) = Node(xpos(a)+xpos(b), ypos(a)+ypos(b), zpos(a)+zpos(b))
Base.:-(a::AbstractNode,b::AbstractNode) = Node(xpos(a)-xpos(b), ypos(a)-ypos(b), zpos(a)-zpos(b))
Base.:*(a::Number,b::AbstractNode) = Node(a*xpos(b), a*ypos(b), a*zpos(b))
Base.:*(b::AbstractNode,a::Number) = Node(a*xpos(b), a*ypos(b), a*zpos(b))
Base.:/(b::AbstractNode,a::Number) = Node(xpos(b)/a, ypos(b)/a, zpos(b)/a)


Base.dot(a::AbstractNode,b::AbstractNode) = muladd(xpos(a), xpos(b), muladd(ypos(a), ypos(b), zpos(a)*zpos(b)))

Base.cross(a::AbstractNode,b::AbstractNode) = Node(ypos(a)*zpos(b) - zpos(a)*ypos(b), zpos(a)*xpos(b) - xpos(a)*zpos(b), xpos(a)*ypos(b) - ypos(a)*xpos(b))

struct Node{T<:Real} <: AbstractNode
  x::T
  y::T
  z::T
end

Node(x,y,z) = Node{Float32}(x,y,z)

#This is kind of a c++ template class. The type of the fields x and y will be the same and a subtype of number.
struct Node2D{T<:Real} <: AbstractNode
  x::T
  y::T
end


#This function returns a zero of the same type of the type contained on the particular Node2D field.
zpos(a::Node2D{T}) where {T} = zero(T)

Base.:+(a::Node2D{T},b::Node2D{T}) where {T} = Node2D{T}(xpos(a)+xpos(b), ypos(a)+ypos(b))
Base.:-(a::Node2D{T},b::Node2D{T}) where {T} = Node2D{T}(xpos(a)-xpos(b), ypos(a)-ypos(b))
Base.:*(a::Number,b::Node2D) = Node2D(a*xpos(b), a*ypos(b))
Base.:*(b::Node2D,a::Number) = Node2D(a*xpos(b), a*ypos(b))
Base.:/(b::Node2D,a::Number) = Node2D(xpos(b)/a, ypos(b)/a)


Base.dot(a::Node2D,b::Node2D) = muladd(xpos(a), xpos(b), ypos(a)*ypos(b))
Base.norm(a::Node2D) = sqrt(muladd(xpos(a), xpos(a), ypos(a)^2))
Base.cross(a::Node2D,b::Node2D) = xpos(a)*ypos(b) - ypos(a)*xpos(b)
