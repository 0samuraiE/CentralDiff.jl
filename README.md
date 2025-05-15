# CentralDiff.jl

*Central difference in Julia*

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://0samuraiE.github.io/CentralDiff.jl/dev/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://0samuraiE.github.io/CentralDiff.jl/dev/)
[![Build Status](https://github.com/0samuraiE/CentralDiff.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/0samuraiE/CentralDiff.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/0samuraiE/CentralDiff.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/0samuraiE/CentralDiff.jl)

CentralDiff.jl is a Julia module for performing central difference on multi-dimensional data. It provides tools to compute numerical derivatives and supports higher-order finite differences.

## Features

- Compute finite differences up to order 20 with zero-cost abstraction.
- Support for multi-dimensional arrays.
- Includes a Simple module for usage that does not require consistency.

## Installation

CentralDiff.jl can be installed with the Julia package manager. From the Julia REPL, type `]` to
enter the Pkg REPL mode and run:
```
pkg> add https://github.com/0samuraiE/CentralDiff.jl.git
```

## Usage
```julia
julia> using CentralDiff

julia> f(x) = x^2;

julia> fc(Order(4), f.(1:4)) # The analytical solution is 2.5^2=6.25
6.25

julia> dfdxc(Order(6), f.(1:6), 1) # The analytical solution is 2*3.5=7.0
7.0

julia> dfdx(Order(4), f.(1:7), 1) # The analytical solution is 2*4.0=8.0
7.999999999999999

julia> d2fdx(Order(4), f.(1:7), 1) # The analytical solution is 2
1.999999999999998

julia> Simple.dfdx(Order(4), f.(1:5), 1) # The analytical solution is 2*3.0=6.0
5.999999999999999

julia> Simple.d2fdx(Order(4), f.(1:5), 1) # The analytical solution is 2
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

julia> fc(Order(2), XAxis(), G, I0) - g̃((x0+dx/2, y0, z0))
-0.0003493436613384304

julia> dfdxc(Order(4), YAxis(), G, I0, dyi) - dg̃dy((x0, y0+dy/2, z0))
-1.4038190432330566e-7

julia> dfdx(Order(6), ZAxis(), G, I0, dzi) - dg̃dz((x0, y0, z0))
-1.7355242798444692e-9

julia> d2fdx(Order(8), ZAxis(), G, I0, dzi) - d2g̃dz((x0, y0, z0))
6.298850330210826e-13

julia> Simple.dfdx(Order(10), XAxis(), G, I0, dxi) - dg̃dx((x0, y0, z0))
-1.1435297153639112e-14

julia> Simple.d2fdx(Order(12), XAxis(), G, I0, dxi) - d2g̃dx((x0, y0, z0))
-1.27675647831893e-15
```
