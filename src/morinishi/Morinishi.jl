# Y. Morinishi, T.S. Lund, O.V. Vasilyev, P. Moin,
# Fully Conservative Higher Order Finite Difference Schemes for Incompressible Flow,
# Journal of Computational Physics, Volume 143, Issue 1, 1998, Pages 90-124,

module Morinishi
export f̄n, δfδxn
using ..CenterDiff

@inline function f̄n(n, F, axis, I)
    Im = I - dI(axis, div(n, 2))
    Ip = I + dI(axis, 1 + div(n, 2))
    (F[Ip] + F[Im]) / 2
end

@inline function δfδxn(n, F, dxi, I)
    Im = I - dI(axis, div(n, 2))
    Ip = I + dI(axis, 1 + div(n, 2))
    (F[Ip] - F[Im]) / n * dxi
end

for n in (1, 3, 5, 7, 9, 11, 13, 15, 17, 19)
    @eval begin
        export $(Symbol(:f̄, n)), $(Symbol(:δfδx, n))
        @inline $(Symbol(:f̄, n))(F, axis, I) = f̄n($n, F, axis, I)
        @inline $(Symbol(:δfδx, n))(F, axis, I) = δfδxn($n, F, axis, I)
    end
end
end
