#ifndef TYPES_HPP
#define TYPES_HPP

template <typename theoryType>
class STT;

/* Type declaration of the ODE systems and the boundary condition minimization function */
template<class theoryType>
struct TYPES{
    using ODE_SYSTEM = void (theoryType::*)(double t, double *y, double *dydt, double data);
    using COST_FUNCTION =  double (STT<theoryType>::*)(double*);    
};
#endif
