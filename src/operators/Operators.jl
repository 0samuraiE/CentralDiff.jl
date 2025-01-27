module Operators
export f̄n, δfδxn
using ..CentralDiff

"""
    f̄n(n, axis, F, I)
    f̄n(n, axis, diff, F, I)

Compute centered average at i+1/2 (default, forward) or i-1/2 (backward) over n points along specified axis.
f̄1, f̄3, f̄5, ..., f̄19 for specific n values available.
"""
@inline function f̄n(n, axis, F, I)
    Im = I - dI(axis, div(n, 2))
    Ip = I + dI(axis, 1 + div(n, 2))
    (F[Ip] + F[Im]) / 2
end

@inline f̄n(n, axis, ::Forward, F, I) = f̄n(n, axis, F, I)
@inline f̄n(n, axis, ::Backward, F, I) = f̄n(n, axis, F, I - dI(axis, 1))

"""
    δfδxn(n, axis, F, I, dxi)
    δfδxn(n, axis, diff, F, I, dxi)

Compute centered difference at i+1/2 (default, forward) or i-1/2 (backward) over n points along specified axis.
δfδx1, δfδx3, δfδx5, ..., δfδx19 for specific n values available.
"""
@inline function δfδxn(n, axis, F, I, dxi)
    Im = I - dI(axis, div(n, 2))
    Ip = I + dI(axis, 1 + div(n, 2))
    (F[Ip] - F[Im]) / n * dxi
end
@inline δfδxn(n, axis, ::Forward, F, I, dxi) = δfδxn(n, axis, F, I, dxi)
@inline δfδxn(n, axis, ::Backward, F, I, dxi) = δfδxn(n, axis, F, I - dI(axis, 1), dxi)

for n in (1, 3, 5, 7, 9, 11, 13, 15, 17, 19)
    @eval begin
        export $(Symbol(:f̄, n)), $(Symbol(:δfδx, n))
        @inline $(Symbol(:f̄, n))(axis, F, I) = f̄n($n, axis, F, I)
        @inline $(Symbol(:f̄, n))(axis, diff, F, I) = f̄n($n, axis, diff, F, I)
        @inline $(Symbol(:δfδx, n))(axis, F, I, dxi) = δfδxn($n, axis, F, I, dxi)
        @inline $(Symbol(:δfδx, n))(axis, diff, F, I, dxi) = δfδxn($n, axis, diff, F, I, dxi)
    end
end
end
