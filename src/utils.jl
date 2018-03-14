abstract type NTupleWrap{N,T} end

Base.@propagate_inbounds Base.getindex(a::NTupleWrap{N,T},i) where {N,T}= Base.getfield(a.ind,i)
Base.eltype(a::Type{NTupleWrap{N,T}}) where {N,T} = T
@inline Base.eltype(a::A) where {A<:NTupleWrap} = eltype(A)
Base.length(a::Type{NTupleWrap{N,T}}) where {N,T} = N
@inline Base.length(a::A) where {A<:NTupleWrap} = length(A)
Base.size(t::NTupleWrap, d) = (d == 1) ? length(t) : throw(ArgumentError("invalid tuple dimension $d"))

Base.start(t::NTupleWrap) =  1
Base.next(t::NTupleWrap,i) = next(t.ind,i)
Base.done(t::NTupleWrap,i) = done(t.ind,i)

Base.:(==)(a::NTupleWrap,b) = a.ind == b
Base.:(==)(b,a::NTupleWrap) = a.ind == b
Base.:(==)(a::NTupleWrap,b::NTupleWrap) = a.ind == b.ind

struct ConstVec{Val} end
Base.getindex(a::Type{ConstVec{Val}}) where {Val} = Val
@inline Base.getindex(a::ConstVec{Val},i) where {Val} = getindex(typeof(a))
Base.eltype(::Type{ConstVec{V}}) where V = typeof(V)
@inline Base.eltype(a::ConstVec) = eltype(typeof(a))

Base.norm(x::Real) = abs(x)