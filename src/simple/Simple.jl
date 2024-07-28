module Simple
export dfdx, d2fdx
using ..CentralDiff
using ..CentralDiff: _generate_terms

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
