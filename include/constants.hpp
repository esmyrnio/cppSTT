#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

namespace CGS
{
    extern const double G;
    extern const double C;
    extern const double MSUN;
    extern const double Length;
    extern const double Time;
    extern const double Density;
};
namespace UNITS
{
    extern const double KAPPA;
    extern const double KSCALE;
};
namespace SURFACE
{
    extern double p_surface;
    extern double e_surface;
    extern const double enthalpy_min;
};
namespace EXECUTION
{
    extern char* call_type;
}
#endif