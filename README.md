# a C++ method to calculate a truncated taylor polynomial for the lorenz system given the parameter values and the start conditions


**Author**: Peter Foelsche |
**Date**: February 2025 |
**Location**: Austin, TX, USA |
**Email**: [peter_foelsche@outlook.com](mailto:peter_foelsche@outlook.com)

## Introduction

[https://en.wikipedia.org/wiki/Lorenz_system](https://en.wikipedia.org/wiki/Lorenz_system)

Good command-line parameters are:

`10 28 2.66666666 0.9 0 0`

This is `sigma`, `rho` and `beta` and `x(t=0)`, `y(t=0)` and `z(t=0)`!

The produced coefficients are polynomials in time for `x(t)`, `y(t)` and `z(t)`.

The number of coefficients can be increased by changing the last template parameter for the `calculate()` function, currently being `3` producing `4` coefficients for each variable. The first coefficient for every variable is obviously identical to the starting conditions passed in via the command-line.

## build

Use the included Visual C++ project file or simply do `g++ -std=c++17 lorenz_taylor.cpp -DNDEBUG -O3`.

## News

I changed the implementation to avoid repeated calculation of the same values referring back to previously calculated time derivatives of the system state.
And I inserted some comments.

## Possible optimizations

I implemented printing out the expressions in order to be able to evaluate them for possible optimizations.
Some easy optimizations I did already implement -- like `something+0==something`.