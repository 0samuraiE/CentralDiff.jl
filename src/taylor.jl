"""
    TaylorMatrix(samples)

Construct a Taylor matrix from the given sample points.
Returns a matrix where element (i,j) is i^j / j!.

# Example
```julia-repl
julia> A = TaylorMatrix(-3//2:3//2)
4×4 Matrix{Rational{Int64}}:
 1  -3//2  9//8  -9//16
 1  -1//2  1//8  -1//48
 1   1//2  1//8   1//48
 1   3//2  9//8   9//16

julia> iA = inv(A)
4×4 Matrix{Rational{Int64}}:
 -1//16   9//16   9//16  -1//16
  1//24  -9//8    9//8   -1//24
  1//2   -1//2   -1//2    1//2
  -1       3      -3       1
```
"""
function TaylorMatrix(samples)
    [i^j / factorial(j) for i in samples, j in 0:length(samples)-1]
end

"""
    _generate_terms(C)

Generate expressions for finite difference or interpolation terms.
Returns a vector of expressions based on coefficient vector C.
"""
function _generate_terms(C)
    N = length(C)
    [:($(Float64(C[n])) * F[I+dI(axis, $(n - cld(N, 2)))]) for n in 1:N]
end

"""
    _shift_and_add(C1, C2)

Compute the convolution of two coefficient vectors C1 and C2.
Returns a vector of length 2N-1, where N is the length of C1 and C2.
"""
function _shift_and_add(C1, C2)
    N = length(C1)
    ret = zeros(eltype(C1), 2N - 1)

    for i in axes(C1, 1), j in axes(C2, 1)
        ret[i+j-1] += C1[i] * C2[j]
    end
    ret
end

"""
    Order{N}()
    Order(N)

Represent the order of a numerical method.
N is an integer constant specifying the order.
"""
struct Order{N} end
Order(N) = Order{N}()

abstract type AbstractDifference end
struct Forward <: AbstractDifference end
struct Backward <: AbstractDifference end

"""
    fc(order, F)
    fc(order, F, axis, I)
    fc(order, F, axis, diff, I)

Compute interpolation of order N along axis at i+1/2 (default, forward) or i-1/2 (backward).
"""
fc(O::Order{N}, F) where {N} = fc(O, F, XAxis(), CartesianIndex(div(N, 2)))
fc(order, F, axis, ::Forward, I) = fc(order, F, axis, I)
fc(order, F, axis, ::Backward, I) = fc(order, F, axis, I - dI(axis, 1))

"""
    dfdxc(order, F, dxi)
    dfdxc(order, F, dxi, axis, I)
    dfdxc(order, F, dxi, axis, diff, I)

Compute first derivative of order N along axis at i+1/2 (default, forward) or i-1/2 (backward).
"""
dfdxc(O::Order{N}, F, dxi) where {N} = dfdxc(O, F, dxi, XAxis(), CartesianIndex(div(N, 2)))
dfdxc(order, F, axis, ::Forward, I) = dfdxc(order, F, dxi, axis, I)
dfdxc(order, F, axis, ::Backward, I) = dfdxc(order, F, dxi, axis, I - dI(axis, 1))

"""
    dfdx(order, F, dxi)
    dfdx(order, F, dxi, axis, I)

Compute first derivative of order N along axis at i.
"""
dfdx(O::Order{N}, F, dxi) where {N} = dfdx(O, F, dxi, XAxis(), CartesianIndex(N))

"""
    d2fdx(order, F, dxi)
    d2fdx(order, F, dxi, axis, I)

Compute second derivative of order N along axis at i.
"""
d2fdx(O::Order{N}, F, dxi) where {N} = d2fdx(O, F, dxi, XAxis(), CartesianIndex(N))

for N in (2, 4, 6, 8, 10, 12, 14, 16, 18, 20)
    samples = Rational{BigInt}.((1-N)//2:(N-1)//2)
    A = TaylorMatrix(samples)
    Ai = inv(A)

    @eval begin
        @inline function fc(::Order{$N}, F, axis, I)
            $(Expr(:call, :+, _generate_terms(@view(Ai[1, :]))...))
        end

        @inline function dfdxc(::Order{$N}, F, dxi, axis, I)
            $(Expr(:call, :+, _generate_terms(@view(Ai[2, :]))...)) * dxi
        end

        @inline function dfdx(::Order{$N}, F, dxi, axis, I)
            $(Expr(:call, :+, _generate_terms(@views(_shift_and_add(Ai[1, :], Ai[2, :])))...)) * dxi
        end

        @inline function d2fdx(::Order{$N}, F, dxi, axis, I)
            $(Expr(:call, :+, _generate_terms(@views(_shift_and_add(Ai[2, :], Ai[2, :])))...)) * dxi^2
        end
    end
end
