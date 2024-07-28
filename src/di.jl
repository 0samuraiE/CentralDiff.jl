"""
    AbstractAxis

An abstract type representing an axis. It serves as the base type for XAxis, YAxis, and ZAxis.
"""
abstract type AbstractAxis end

struct XAxis <: AbstractAxis end
struct YAxis <: AbstractAxis end
struct ZAxis <: AbstractAxis end

"""
    dI{T<:AbstractAxis}

A structure representing movement along an axis.

Fields:
- `axis::T`: The axis along which movement occurs
- `stride::Int`: The amount of movement (stride)

`T` must be a subtype of `AbstractAxis`.

# Examples
```jldoctest
julia> CartesianIndex(1, 1) + dI(XAxis(), 2)
CartesianIndex(3, 1)

julia> CartesianIndex(1, 1) + dI(YAxis(), 2)
CartesianIndex(1, 3)

julia> CartesianIndex(1, 1, 1) + dI(ZAxis(), 3)
CartesianIndex(1, 1, 4)
"""
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
-(x::CartesianIndex, y::dI) = +(x, dI(y.axis, -y.stride))
-(x::dI, y::CartesianIndex) = +(y, x)
