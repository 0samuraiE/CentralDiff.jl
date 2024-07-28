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

"""
    countorder(f, hs...)

Compute orders of convergence for a given function across multiple grid spacings.

# Arguments
- `f`: A function that computes error(s) for a given grid spacing.
- `hs...`: A variadic list of grid spacings.

# Returns
Either a vector of computed orders (for single error measure) or a NamedTuple of 
computed orders (for multiple error measures).

# Description
This function is a high-level interface for computing orders of convergence. It applies 
the given function `f` to each grid spacing in `hs`, then computes the orders of convergence 
based on the returned errors. It can handle both single and multiple error measures.
"""
function countorder(f, hs...)
    errslist = map(f, hs)
    _countorder(hs, errslist)
end