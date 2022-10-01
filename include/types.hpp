#ifndef TYPES_HPP
#define TYPES_HPP

// forward declaration
template <typename theoryType>
class STT;
template<class theoryType>
struct TYPES{
    using ODE_SYSTEM = void (theoryType::*)(double t, double *y, double *dydt, double data);
    using COST_FUNCTION =  double (STT<theoryType>::*)(double*);    
};
#endif