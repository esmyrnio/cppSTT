#include "../include/constants.hpp"
#include <math.h>

char* EXECUTION::call_type = "get_MR";

const double CGS::G = 6.6732e-8;
const double CGS::C = 2.9979e+10;
const double CGS::MSUN = 1.989e+33;
const double CGS::Length = G*MSUN/pow(C,2);
const double CGS::Density = MSUN/pow(Length,3);
const double UNITS::KSCALE = 1.112668301525780e-36;
const double UNITS::KAPPA = 1.346790806509621e+13;
const double SURFACE::enthalpy_min = 1.0/(CGS::C*CGS::C);