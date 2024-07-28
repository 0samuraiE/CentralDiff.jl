"""
    countorder(f, hs...)

Compute convergence orders for a function over multiple grid spacings.
Returns computed orders as a NamedTuple or vector.
"""
function countorder(f, hs...)
    errslist = map(f, hs)
    _countorder(hs, errslist)
end

function _compute_orders(hs, errs)
    log.(hs[2:end] ./ hs[1:end-1], errs[2:end] ./ errs[1:end-1])
end

function _countorder(hs, errslist::NTuple{N,NamedTuple}) where {N}
    keys = propertynames(errslist[begin])
    orderlist = map(keys) do key
        es = getfield.(errslist, key)
        _compute_orders(hs, es)
    end
    NamedTuple{keys}(orderlist)
end

function _countorder(hs, errslist::NTuple{N,Float64}) where {N}
    es = errslist
    _compute_orders(hs, es)
end
