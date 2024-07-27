module CenterDiff
export Axis, XAxis, YAxis, ZAxis, dI, Order
export TaylorMatrix, f̄c, dfdxc, dfdx, d2fdx, Simple, countorder

import Base: +

struct Order{N} end
Order(N) = Order{N}()

abstract type AbstractAxis end
struct XAxis <: AbstractAxis end
struct YAxis <: AbstractAxis end
struct ZAxis <: AbstractAxis end

struct dI{T<:AbstractAxis}
    axis::T
    stride::Int
end

+(x::CartesianIndex{1}, y::dI{XAxis}) = x + CartesianIndex(y.stride)
+(x::CartesianIndex{2}, y::dI{XAxis}) = x + CartesianIndex(y.stride, 0)
+(x::CartesianIndex{2}, y::dI{YAxis}) = x + CartesianIndex(0, y.stride)
+(x::CartesianIndex{3}, y::dI{XAxis}) = x + CartesianIndex(y.stride, 0, 0)
+(x::CartesianIndex{3}, y::dI{YAxis}) = x + CartesianIndex(0, y.stride, 0)
+(x::CartesianIndex{3}, y::dI{ZAxis}) = x + CartesianIndex(0, 0, y.stride)
+(x::dI, y::CartesianIndex) = +(y, x)

function TaylorMatrix(samples)
    [i^j / factorial(j) for i in samples, j in 0:length(samples)-1]
end

function _generate_terms(C)
    N = length(C)
    [:($(Float64(C[n])) * F[I+dI(axis, $(n - cld(N, 2)))]) for n in 1:N]
end

function _shift_and_add(C1, C2)
    N = length(C1)
    ret = zeros(eltype(C1), 2N - 1)

    for i in axes(C1, 1), j in axes(C2, 1)
        ret[i+j-1] += C1[i] * C2[j]
    end
    ret
end

f̄c(O::Order{N}, F) where {N} = f̄c(O, XAxis(), F, CartesianIndex(div(N, 2)))
dfdxc(O::Order{N}, F, dxi) where {N} = dfdxc(O, XAxis(), F, CartesianIndex(div(N, 2)), dxi)
dfdx(O::Order{N}, F, dxi) where {N} = dfdx(O, XAxis(), F, CartesianIndex(N), dxi)
d2fdx(O::Order{N}, F, dxi) where {N} = d2fdx(O, XAxis(), F, CartesianIndex(N), dxi)

for N in (2, 4, 6, 8, 10, 12, 14, 16, 18, 20)
    samples = Rational{BigInt}.((1-N)//2:(N-1)//2)
    A = TaylorMatrix(samples)
    Ai = inv(A)

    @eval begin
        @inline function f̄c(::Order{$N}, axis, F, I)
            $(Expr(:call, :+, _generate_terms(@view(Ai[1, :]))...))
        end

        @inline function dfdxc(::Order{$N}, axis, F, I, dxi)
            $(Expr(:call, :+, _generate_terms(@view(Ai[2, :]))...)) * dxi
        end

        @inline function dfdx(::Order{$N}, axis, F, I, dxi)
            $(Expr(:call, :+, _generate_terms(@views(_shift_and_add(Ai[1, :], Ai[2, :])))...)) * dxi
        end

        @inline function d2fdx(::Order{$N}, axis, F, I, dxi)
            $(Expr(:call, :+, _generate_terms(@views(_shift_and_add(Ai[2, :], Ai[2, :])))...)) * dxi^2
        end
    end
end

module Simple
    using ..CenterDiff
    using ..CenterDiff: _generate_terms

    @inline dfdx(O::Order{N}, F, dxi) where {N} = dfdx(O, XAxis(), F, CartesianIndex(div(N, 2) + 1), dxi)
    @inline d2fdx(O::Order{N}, F, dxi) where {N} = d2fdx(O, XAxis(), F, CartesianIndex(div(N, 2) + 1), dxi)

    for N in (2, 4, 6, 8, 10, 12, 14, 16, 18, 20)
        samples = Rational{BigInt}.(-div(N, 2):div(N, 2))
        A = TaylorMatrix(samples)
        Ai = inv(A)

        @eval begin
            @inline function dfdx(::Order{$N}, axis, F, I, dxi)
                $(Expr(:call, :+, _generate_terms(@view(Ai[2, :]))...)) * dxi
            end

            @inline function d2fdx(::Order{$N}, axis, F, I, dxi)
                $(Expr(:call, :+, _generate_terms(@view(Ai[3, :]))...)) * dxi^2
            end
        end
    end
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

function countorder(f, hs...)
    errslist = map(f, hs)
    _countorder(hs, errslist)
end
end
