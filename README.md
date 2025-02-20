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

The number of coefficients can be changed by changing the constant `MAX_ORDER`. The first coefficient for every variable (`X`, `Y`, `Z`) is obviously identical to the starting conditions passed in via the command-line.

The output for the above commandline is now the following:
```
peter@M4700:~/lorenz_taylor$ ./lorentz_taylor.exe 10 28 2.66666666 0.9 0 0
X<1>=(sigma*(Y<0>-X<0>))
Y<1>=((X<0>*(rho-Z<0>))-Y<0>)
Z<1>=((X<0>*Y<0>)-(beta*Z<0>))
X<2>=(sigma*(Y<1>-X<1>))
Y<2>=(((X<1>*(rho-Z<0>))-(X<0>*Z<1>))-Y<1>)
Z<2>=(((X<1>*Y<0>)+(X<0>*Y<1>))-(beta*Z<1>))
X<3>=(sigma*(Y<2>-X<2>))
Y<3>=((((X<2>*(rho-Z<0>))-(X<1>*Z<1>))-((X<1>*Z<1>)+(X<0>*Z<2>)))-Y<2>)
Z<3>=((((X<2>*Y<0>)+(X<1>*Y<1>))+((X<1>*Y<1>)+(X<0>*Y<2>)))-(beta*Z<2>))
(0.9, d/dX0=1)
(-9, d/dX0=-10, d/dX1=10)
(171, d/dX0=190, d/dX1=-55, d/dX2=-4.5)
(-1032, d/dX0=-1146.67, d/dX1=650.317, d/dX2=35.5)

(0, d/dX1=1)
(25.2, d/dX0=28, d/dX1=-1, d/dX2=-0.9)
(-138.6, d/dX0=-154, d/dX1=140.095, d/dX2=6.15)
(1638.8, d/dX0=1813.33, d/dX1=-555.487, d/dX2=-109.995)

(0, d/dX2=1)
(0, d/dX1=0.9, d/dX2=-2.66667)
(11.34, d/dX0=25.2, d/dX1=-6.15, d/dX2=3.15056)
(-127.26, d/dX0=-282.8, d/dX1=191.495, d/dX2=1.74451)

peter@M4700:~/lorenz_taylor$
```

The coefficients are 3 sets for `x`, `y`, `z` and they are 4 values each (`MAX_ORDER` + 1). The value is the first value after the opening parenthesis and the labled values following are derivatives wrt the initial conditions specified.
## build

Use the included Visual C++ project file or simply do `g++ -std=c++17 lorenz_taylor.cpp -DNDEBUG -O3 -I $BOOST_ROOT/include`.
I also added a `Makefile`.
If you do not use either the `Makefile` nor the Visual C++ solution file, you will have to initialize the submodule by hand:
```
git submodule init ctaylor
git submodule update --init --recursive ctaylor
```

## News

Now this code references my github repository [excessphase/ctaylor](https://github.com/ExcessPhase/ctaylor) for calculating the derivative of the taylor coefficients wrt the initial conditions of the state variables.


I changed the implementation to avoid repeated calculation of the same values referring back to previously calculated time derivatives of the system state.
And I inserted some comments.

## Possible optimizations

I implemented printing out the expressions in order to be able to evaluate them for possible optimizations.
Some easy optimizations I did already implement -- like `something+0==something`.

## Looking for a job as a C++ Software Engineer

BTW -- if you know some company which is trying to hire a C++ Software Engineer: I'm on the market since some months already!
