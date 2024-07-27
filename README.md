# CenterDiff.jl

*Central difference in Julia*

[![Build Status](https://github.com/0samuraiE/CenterDiff.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/0samuraiE/CenterDiff.jl/actions/workflows/CI.yml?query=branch%3Amaster)

CenterDiff.jl is a Julia module for performing central difference on multi-dimensional data. It provides tools to compute numerical derivatives and supports higher-order finite differences.

## Features

- Compute finite differences up to order 20.
- Support for multi-dimensional arrays.
- Includes a Simple module for usage that does not require consistency.

## Installation

CenterDiff.jl can not be installed with the Julia package manager. From the Julia REPL, type `]` to
enter the Pkg REPL mode and run:
```
pkg> add git@github.com:0samuraiE/CenterDiff.jl.git
```
