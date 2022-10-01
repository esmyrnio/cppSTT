#!/bin/bash

swig -python -c++ pySTT.i
g++ -std=c++11 -c -O3 -fpic ../src/main.cpp ../src/constants.cpp ../src/eos.cpp
# here use path to your python source
g++ -std=c++11 -c -O3 -fpic pySTT_wrap.cxx -I/usr/include/python2.7
g++ -shared -O3 pySTT_wrap.o main.o constants.o eos.o -o _pySTT.so
