# CenterDiff.jl

*Central difference in Julia*

[![Build Status](https://github.com/0samuraiE/CenterDiff.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/0samuraiE/CenterDiff.jl/actions/workflows/CI.yml?query=branch%3Amaster)

CenterDiff.jl is a Julia module for performing central difference on multi-dimensional data. It provides tools to compute numerical derivatives and supports higher-order finite differences.

## Features

- Compute finite differences up to order 20 with zero-cost abstraction.
- Support for multi-dimensional arrays.
- Includes a Simple module for usage that does not require consistency.

## Installation

CenterDiff.jl can not be installed with the Julia package manager. From the Julia REPL, type `]` to
enter the Pkg REPL mode and run:
```
pkg> add https://github.com/0samuraiE/CenterDiff.jl.git
```

## Usage
```julia
julia> using CenterDiff

julia> f(x) = x^2;

julia> f̄c(Order(4), f.(1:4)) # The analytical solution is 2.5^2=6.25
6.25
julia> f̄c(Order(4), f.(1:4), XAxis(), CartesianIndex(div(4, 2))) # equivalent
6.25

julia> dfdxc(Order(6), f.(1:6), 1) # The analytical solution is 2*3.5=7.0
7.0
julia> dfdxc(Order(6), [1, 2, 3, 4, 5, 6], 1, XAxis(), CartesianIndex(div(6, 2))) # equivalent
7.0

julia> dfdx(Order(4), f.(1:7), 1) # The analytical solution is 2*4.0=8.0
7.999999999999999
julia> dfdx(Order(4), f.(1:7), 1, XAxis(), CartesianIndex(4)) # equivalent
7.999999999999999

julia> d2fdx(Order(4), f.(1:7), 1) # The analytical solution is 2
1.999999999999998
julia> d2fdx(Order(4), f.(1:7), 1, XAxis(), CartesianIndex(4)) # equivalent
1.999999999999998

julia> Simple.dfdx(Order(4), f.(1:5), 1) # The analytical solution is 2*3.0=6.0
5.999999999999999
julia> Simple.dfdx(Order(4), f.(1:5), 1, XAxis(), CartesianIndex(3)) # equivalent
5.999999999999999

julia> Simple.d2fdx(Order(4), f.(1:5), 1) # The analytical solution is 2
1.9999999999999991
julia> Simple.d2fdx(Order(4), f.(1:5), 1, XAxis(), CartesianIndex(3)) # equivalent
1.9999999999999991
```

## Advanced Usage
```julia
julia> g̃((x, y, z)) = sin(x) * sin(y) * sin(z);
       dg̃dx((x, y, z)) = cos(x) * sin(y) * sin(z);
       dg̃dy((x, y, z)) = sin(x) * cos(y) * sin(z);
       dg̃dz((x, y, z)) = sin(x) * sin(y) * cos(z);
       d2g̃dx((x, y, z)) = -sin(x) * sin(y) * sin(z);
       d2g̃dy((x, y, z)) = -sin(x) * sin(y) * sin(z);
       d2g̃dz((x, y, z)) = -sin(x) * sin(y) * sin(z);
       X = Y = Z = range(0, 2π; length=64);
       dx = dy = dz = step(X); dxi = dyi = dzi = 1 / dx;
       MESH = collect(Iterators.product(X, Y, Z));
       G = g̃.(MESH);
       I0 = CartesianIndex(8, 8, 8);
       x0, y0, z0 = MESH[I0];

julia> f̄c(Order(2), G, XAxis(), I0) - g̃((x0+dx/2, y0, z0))
-0.0003493436613384304

julia> dfdxc(Order(4), G, dyi, YAxis(), I0) - dg̃dy((x0, y0+dy/2, z0))
-1.4038190432330566e-7

julia> dfdx(Order(6), G, dzi, ZAxis(), I0) - dg̃dy((x0, y0, z0))
-1.7355242798444692e-9

julia> d2fdx(Order(8), G, dzi, ZAxis(), I0) - d2g̃dz((x0, y0, z0))
6.298850330210826e-13

julia> Simple.dfdx(Order(10), G, dxi, XAxis(), I0) - dg̃dy((x0, y0, z0))
-1.1435297153639112e-14

julia> Simple.d2fdx(Order(12), G, dxi, XAxis(), I0) - d2g̃dz((x0, y0, z0))
-1.27675647831893e-15
```
