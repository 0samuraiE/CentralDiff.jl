# Y. Morinishi, T.S. Lund, O.V. Vasilyev, P. Moin,
# Fully Conservative Higher Order Finite Difference Schemes for Incompressible Flow,
# Journal of Computational Physics, Volume 143, Issue 1, 1998, Pages 90-124,

module Morinishi
export f̄n, δfδxn
using ..CenterDiff

"""
    f̄n(n, F, axis, I)
    f̄n(n, F, axis, diff, I)

Compute centered average at i+1/2 (default, forward) or i-1/2 (backward) over n points along specified axis.
f̄1, f̄3, f̄5, ..., f̄19 for specific n values available.
"""
@inline function f̄n(n, F, axis, I)
    Im = I - dI(axis, div(n, 2))
    Ip = I + dI(axis, 1 + div(n, 2))
    (F[Ip] + F[Im]) / 2
end

@inline f̄n(n, F, axis, ::Forward, I) = f̄n(n, F, axis, I)
@inline f̄n(n, F, axis, ::Backward, I) = f̄n(n, F, axis, I - dI(axis, 1))

"""
    δfδxn(n, F, dxi, I)
    δfδxn(n, F, dxi, axis, diff, I)

Compute centered difference at i+1/2 (default, forward) or i-1/2 (backward) over n points along specified axis.
δfδx1, δfδx3, δfδx5, ..., δfδx19 for specific n values available.
"""
@inline function δfδxn(n, F, dxi, axis, I)
    Im = I - dI(axis, div(n, 2))
    Ip = I + dI(axis, 1 + div(n, 2))
    (F[Ip] - F[Im]) / n * dxi
end
@inline δfδxn(n, F, dxi, axis, ::Forward, I) = δfδxn(n, F, dxi, axis, I)
@inline δfδxn(n, F, dxi, axis, ::Backward, I) = δfδxn(n, F, dxi, axis, I - dI(axis, 1))

for n in (1, 3, 5, 7, 9, 11, 13, 15, 17, 19)
    @eval begin
        export $(Symbol(:f̄, n)), $(Symbol(:δfδx, n))
        @inline $(Symbol(:f̄, n))(F, axis, I) = f̄n($n, F, axis, I)
        @inline $(Symbol(:f̄, n))(F, axis, diff, I) = f̄n($n, F, axis, diff, I)
        @inline $(Symbol(:δfδx, n))(F, dxi, axis, I) = δfδxn($n, F, dxi, axis, I)
        @inline $(Symbol(:δfδx, n))(F, dxi, axis, diff, I) = δfδxn($n, F, dxi, axis, diff, I)
    end
end
end
