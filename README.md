# cppSTT

## Description

Program for obtaining neutron star models for tabulated equations of state[^1], in the context of Damour-Esposito-Farese (DEF) scalar tensor theory model and R^2 gravity model. The system of ODEs is solved by implementing a C++11 version of LSODA[^2]. Additionally, a C++ version of Nelder-Mead method[^3] is used for the optimization of the initial conditions as well as some modules from the BOOST C++ library[^4]. The code prints the model's parameters e.g. mass, radius, scalar charge, scalar at center and boundary minimization accuracy, for the corresponding input parameters.

The C++ code is wrapped using SWIG[^5], into a python library which can be imported and used in any way, as shown in swig/test.py.

## Usage

Compile, using the provided makefile. Type "make" in the main directory, where the makefile is also located, to create the executionable file.

STT.exe takes 5 inputs in-turn. The parameters are specified using the following flags:

1. **-m** *theory* (DEF or R2)
2. **-f** *eos name* (EoS file name).
3. **-c** *coupling* (model's coupling constant).
4. **-e** *central density* (central energy density in CGS/10^15).

For example,

```
./STT -m DEF -f ppsly4.cold.rns1.1.txt -c -5.0 -e 1.5
```
```
./STT -m R2 -f ppsly4.cold.rns1.1.txt -c 20.0 -e 1.5
```
*typical input parameter space:*
   - *central density ~ 0.5-4.0 (10^15 gr/cm^3)*
   - *coupling ~ (-5.0,-4.5) (DEF), (5,10^3) (R^2)*

GR solutions may also be obtained simply by choosing "DEF" as theory and setting the coupling to zero.

[^1]:https://ui.adsabs.harvard.edu/abs/2021PhRvD.103l3004B/abstract (table IX)
[^2]:https://github.com/dilawar/libsoda-cxx
[^3]:https://people.sc.fsu.edu/~jburkardt/cpp_src/asa047/asa047.html
[^4]:https://github.com/boostorg/boost
[^5]:https://github.com/swig/swig
