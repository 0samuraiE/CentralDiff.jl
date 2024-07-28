module Simple
export dfdx, d2fdx
using ..CenterDiff
using ..CenterDiff: _generate_terms

@inline dfdx(O::Order{N}, F, dxi) where {N} = dfdx(O, F, dxi, XAxis(), CartesianIndex(div(N, 2) + 1))
@inline d2fdx(O::Order{N}, F, dxi) where {N} = d2fdx(O, F, dxi, XAxis(), CartesianIndex(div(N, 2) + 1))

for N in (2, 4, 6, 8, 10, 12, 14, 16, 18, 20)
    samples = Rational{BigInt}.(-div(N, 2):div(N, 2))
    A = TaylorMatrix(samples)
    Ai = inv(A)

    @eval begin
        @inline function dfdx(::Order{$N}, F, dxi, axis, I)
            $(Expr(:call, :+, _generate_terms(@view(Ai[2, :]))...)) * dxi
        end

        @inline function d2fdx(::Order{$N}, F, dxi, axis, I)
            $(Expr(:call, :+, _generate_terms(@view(Ai[3, :]))...)) * dxi^2
        end
    end
end
end
