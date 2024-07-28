"""
    TaylorMatrix(samples)

Construct a Taylor matrix from the given samples.

# Arguments
- `samples`: A vector of sample points.

# Returns
A matrix where each element (i,j) is calculated as i^j / j!, 
where i is a sample point and j ranges from 0 to (length(samples) - 1).

# Description
This function creates a matrix based on the Taylor series expansion. 
Each row corresponds to a sample point, and each column represents 
a term in the Taylor series up to the (n-1)th derivative, where n is 
the number of samples.

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
Generate a list of expressions representing terms for computing approximations or derivatives.

# Arguments
- `C`: A vector of coefficients.

# Returns
A vector of expressions, each representing a term.

# Details
This function creates expressions for each term. 
It uses the coefficients in `C` to create terms that will later be summed.

Each term is of the form: `coefficient * F[I + dI(axis, offset)]`
where:
- `coefficient` is an element from `C`
- `F` is assumed to be the array being computed
- `I` is assumed to be the current index
- `axis` is assumed to be defined in the context where this function is called
- `offset` is calculated based on the position of the coefficient in `C`

# Examples
```julia-repl
julia> C = [1.0, 2.0, 1.0];

julia> terms = _generate_terms(C)
3-element Vector{Expr}:
 :(1.0 * F[I + dI(axis, -1)])
 :(2.0 * F[I + dI(axis, 0)])
 :(1.0 * F[I + dI(axis, 1)])
```
"""
function _generate_terms(C)
    N = length(C)
    [:($(Float64(C[n])) * F[I+dI(axis, $(n - cld(N, 2)))]) for n in 1:N]
end

"""
    _shift_and_add(C1, C2)

Perform a convolution-like operation on two coefficient vectors.

# Arguments
- `C1`: First vector of coefficients.
- `C2`: Second vector of coefficients.

# Returns
A vector containing the result of the shift-and-add operation.

# Description
This function performs a discrete convolution-like operation on two input vectors. 
It is particularly useful for combining finite difference coefficients or stencils.
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
    Order{N}

A type representing the order of a numerical method or approximation.

# Type Parameters
- `N`: An integer constant representing the order.
"""
struct Order{N} end
Order(N) = Order{N}()

"""
    fc(O::Order{N}, F) where {N}
    fc(O::Order{N}, F, axis::AbstractAxis, I::CartesianIndex) where {N}

Compute the centered interpolation at i+1/2 for the given order and data.

# Arguments
- `O::Order{N}`: Order of the interpolation.
- `F`: Data array to interpolate.
- `axis::AbstractAxis`: Axis along which to perform the interpolation (default: `XAxis()`).
- `I::CartesianIndex`: Base index for the operation (default: `CartesianIndex(div(N, 2))`).

# Returns
The interpolated value at the i+1/2 position.

# Description
This function computes a centered interpolation at the midpoint between grid points (i+1/2) 
for the given order `N`. The 'c' in the function name indicates that it's a centered operation at i+1/2.
"""
fc(O::Order{N}, F) where {N} = fc(O, F, XAxis(), CartesianIndex(div(N, 2)))

"""
    dfdxc(O::Order{N}, F, dxi) where {N}
    dfdxc(O::Order{N}, F, dxi, axis::AbstractAxis, I::CartesianIndex) where {N}

Compute the centered first derivative at i+1/2 for the given order and data.

# Arguments
- `O::Order{N}`: Order of the derivative approximation.
- `F`: Data array to differentiate.
- `dxi`: Grid spacing or the reciprocal of Δx.
- `axis::AbstractAxis`: Axis along which to perform the differentiation (default: `XAxis()`).
- `I::CartesianIndex`: Base index for the operation (default: `CartesianIndex(div(N, 2))`).

# Returns
The approximated first derivative at the i+1/2 position.

# Description
This function computes a centered finite difference approximation of the first derivative 
at the midpoint between grid points (i+1/2) for the given order `N`. The 'c' in the 
function name indicates that it's a centered operation at i+1/2.
"""
dfdxc(O::Order{N}, F, dxi) where {N} = dfdxc(O, F, dxi, XAxis(), CartesianIndex(div(N, 2)))

"""
    dfdx(O::Order{N}, F, dxi) where {N}
    dfdx(O::Order{N}, F, dxi, axis::AbstractAxis, I::CartesianIndex) where {N}

Compute the centered first derivative at i for the given order and data.

# Arguments
- `O::Order{N}`: Order of the derivative approximation.
- `F`: Data array to differentiate.
- `dxi`: Grid spacing or the reciprocal of Δx.
- `axis::AbstractAxis`: Axis along which to perform the differentiation (default: `XAxis()`).
- `I::CartesianIndex`: Base index for the operation (default: `CartesianIndex(N)`).

# Returns
The approximated first derivative at the i position.

# Description
This function computes a centered finite difference approximation of the first derivative 
at the grid points (i) for the given order `N`.
"""
dfdx(O::Order{N}, F, dxi) where {N} = dfdx(O, F, dxi, XAxis(), CartesianIndex(N))

"""
    d2fdx(O::Order{N}, F, dxi) where {N}
    d2fdx(O::Order{N}, F, dxi, axis::AbstractAxis, I::CartesianIndex) where {N}

Compute the centered second derivative at i for the given order and data.

# Arguments
- `O::Order{N}`: Order of the derivative approximation.
- `F`: Data array to differentiate.
- `dxi`: Grid spacing or the reciprocal of Δx.
- `axis::AbstractAxis`: Axis along which to perform the differentiation (default: `XAxis()`).
- `I::CartesianIndex`: Base index for the operation (default: `CartesianIndex(N)`).

# Returns
The approximated second derivative at the i position.

# Description
This function computes a centered finite difference approximation of the second derivative 
at the grid points (i) for the given order `N`.
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
