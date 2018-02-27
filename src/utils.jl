abstract type NTupleWrap{N,T} end

Base.@propagate_inbounds Base.getindex(a::NTupleWrap{N,T},i) where {N,T}= Base.getfield(a.ind,i)
Base.eltype(a::NTupleWrap{N,T}) where {N,T} = T
Base.length(a::NTupleWrap{N,T}) where {N,T} = N
Base.size(t::NTupleWrap, d) = (d == 1) ? length(t) : throw(ArgumentError("invalid tuple dimension $d"))

Base.start(t::NTupleWrap) =  1
Base.next(t::NTupleWrap,i) = next(t.ind,i)
Base.done(t::NTupleWrap,i) = done(t.ind,i)

Base.:(==)(a::NTupleWrap,b) = a.ind == b
Base.:(==)(b,a::NTupleWrap) = a.ind == b
Base.:(==)(a::NTupleWrap,b::NTupleWrap) = a.ind == b.ind

struct ConstVec{Val} end
Base.getindex(a::ConstVec{Val},i) where {Val} = Val
