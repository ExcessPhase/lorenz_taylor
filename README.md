# a method to calculate a truncated taylor polynomial for the lorenz system given the parameter values and the start condiions


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

The number of coefficients can be increased by changing the last template parameter.for the `call()` function, currently being `3` producing `4` coefficients for each variable.