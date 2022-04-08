# STT_solvers

## Description

Constructing neutron star models in the context of Damour-Esposito-Farese (DEF) scalar-tensor model and R^2 - gravity model. The system of ODEs is solved by implementing the c++11 version of LSODA (https://github.com/dilawar/libsoda-cxx). Additionally, a c++ version of Nelder-Mead method (https://people.sc.fsu.edu/~jburkardt/cpp_src/asa047/asa047.html) is used for the optimization of the initial conditions. 

## Usage

Compile, using the provided makefile. Simply type "make" in the main directory, where the makefile is also located.

run ./DEF or ./R2

DEF.exe and R2.exe take 4 inputs in-turn:

1. The EoS file name (e.g. eosSLY.txt).
2. The coupling constant (e.g. -5.0 for DEF and 100.0 for R^2).
3. The central energy density in CGS/10^15 (e.g 1.2).
4. The print option 0 or 1:
    -  0: Prints minimization info along with gravitational mass M, radius R and scalar charge.
    -  1: Prints (0) along with the distance, metric, scalar, energy density and pressure profiles.
