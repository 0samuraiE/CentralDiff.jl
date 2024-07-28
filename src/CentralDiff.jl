module CentralDiff
export Axis, XAxis, YAxis, ZAxis, dI, Order
export TaylorMatrix, fc, dfdxc, dfdx, d2fdx, Simple, Morinishi

include("di.jl")
include("taylor.jl")
include("morinishi/Morinishi.jl")
include("simple/Simple.jl")
include("uitl.jl")
end