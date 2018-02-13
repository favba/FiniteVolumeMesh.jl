abstract type NTupleWrap{N,T} end

Base.@propagate_inbounds Base.getindex(a::NTupleWrap,i) = getindex(a.ind,i)
Base.eltype(a::NTupleWrap{N,T}) where {N,T} = T
Base.length(a::NTupleWrap{N,T}) where {N,T} = N
Base.size(t::NTupleWrap, d) = (d == 1) ? length(t) : throw(ArgumentError("invalid tuple dimension $d"))

Base.:(==)(a::NTupleWrap,b) = a.ind == b
Base.:(==)(b,a::NTupleWrap) = a.ind == b
Base.:(==)(a::NTupleWrap,b::NTupleWrap) = a.ind == b.ind