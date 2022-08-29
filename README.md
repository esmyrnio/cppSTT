# STT_solvers

## Description

Program for obtaining neutron star models in the context of Damour-Esposito-Farese (DEF) scalar-tensor theory model and R^2 - gravity model. The system of ODEs is solved by implementing a C++11 version of LSODA^[1]. Additionally, a C++ version of Nelder-Mead method ^[2] is used for the optimization of the initial conditions as well as some modules from the BOOST C++ library ^[3]. 

The C++ code is wrapped using SWIG[^4], into a python library which can be imported and used in any way, as shown in swig/test.py.

## Usage

Compile, using the provided makefile. Type "make" in the main directory, where the makefile is also located, to create the executionable file.

STT.exe takes 5 inputs in-turn. The parameters are specified using the following flags:

1. **-m** *theory* (DEF or R2)
2. **-f** *eos name* (EoS file name).
3. **-c** *coupling* (model's coupling constant).
4. **-e** *central density* (central energy density in CGS/10^15).

For example,

```
./STT -m DEF -f ppsly4.cold.rns1.1.txt -c -5.0 -e 1.5 -p 0
```
```
./STT -m R2 -f ppsly4.cold.rns1.1.txt -c 20.0 -e 1.5 -p 0
```

The code prints the model's parameters e.g. mass, radius, scalar charge, scalar at center and boundary minimization accuracy, for the corresponding input parameters.

[^1]:https://github.com/dilawar/libsoda-cxx
[^2]:https://people.sc.fsu.edu/~jburkardt/cpp_src/asa047/asa047.html
[^3]:https://github.com/boostorg/boost
[^4]:https://github.com/swig/swig
