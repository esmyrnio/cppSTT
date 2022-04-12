# STT_solvers

## Description

Constructing neutron star models in the context of Damour-Esposito-Farese (DEF) scalar-tensor model and R^2 - gravity model. The system of ODEs is solved by implementing a C++11 version of LSODA (https://github.com/dilawar/libsoda-cxx). Additionally, a C++ version of Nelder-Mead method (https://people.sc.fsu.edu/~jburkardt/cpp_src/asa047/asa047.html) is used for the optimization of the initial conditions. 

## Usage

Compile, using the provided makefile. Simply type "make" in the main directory, where the makefile is also located, to create the executionables.

DEF.exe and R2.exe take 4 inputs in-turn. The parameters are specified using the following flags:

1. **-f** *eos_name* (The EoS file name).
2. **-c** *coupling* (The EGB coupling constant in km^2).
3. **-e** *central_density* (The central energy density in CGS/10^15).
4. **-p** *The print option 0 or 1*:
    -  0: Prints gravitational mass M and radius R.
    -  1: Prints (0) along with the distance, metric, scalar, energy density and pressure profiles.

For example,

*./DEF -f eosSLY.txt -c -5.0 -e 1.5 -p 0*

*./R2 -f eosSLY.txt -c 10.0 -e 1.5 -p 0*
