#ifndef THEORY_HPP
#define THEORY_HPP

#include <vector>
#include <string>
#include "eos.hpp"
#include "constants.hpp"
#include <math.h>

using Vector = std::vector<double>;
using Array = double*;

struct DEF
{
    DEF(std::string,double foo,double);
    ~DEF();
    double A_of_phi(double, double);
    double a_of_phi(double,double);
    double V_of_phi(double,double);
    double dV_of_phi(double, double);
    double dnudr_inf(double,double,double,double);
    void system(double,Array,Array,double);
    void computeMass(double&,double,double,double, double);
    void computeRadius(double&, double,double,double);
    void computeScalarCharge(double&,double);
    double scalar_start(double eps_c);
    double r_max,r_max_surface,RDIV;
    double atol,rtol,reqmin;
    double A, alpha,eps;
    double startingPointX;
    double startingPointY;
    double stepX;
    double stepY;
    EOS eos;
};

struct R2
{
    R2(std::string,double,double foo);
    ~R2();
    double A_of_phi(double, double);
    double a_of_phi(double,double);
    double V_of_phi(double,double);
    double dV_of_phi(double, double);
    double dnudr_inf(double,double,double,double);
    void system(double,Array,Array,double);
    void computeMass(double&,double,double,double, double);
    void computeRadius(double&, double,double,double);
    void computeScalarCharge(double&,double);
    double set_rmax(double);
    double r_max,r_max_surface,RDIV;
    double atol,rtol,reqmin;
    double A, V,dV,alpha,eps;
    double startingPointX;
    double startingPointY;
    double stepX;
    double stepY;
    EOS eos;
};
#endif

DEF::DEF(std::string eos_name, double foo, double eps_c):atol(1e-15),rtol(1e-06),reqmin(1e-15),
r_max(1e5),r_max_surface(15.0),RDIV(5001),startingPointX(-0.8),startingPointY(scalar_start(eps_c)),
stepX(0.2),stepY(-0.01), eos(eos_name){};

DEF::~DEF(){};

double DEF::scalar_start(double eps_c)
{
    if(eps_c>2.1)
    {
        return 0.0;
    }
    else
    {
        return -0.1;
    }
    
}
double DEF::a_of_phi(double phi, double beta)
{
    return beta*phi;
}
double DEF::V_of_phi(double phi, double beta)
{
    return 0;
}
double DEF::dV_of_phi(double phi, double beta)
{
    return 0;
}
double DEF::A_of_phi(double phi, double beta)
{
    return exp(0.5*beta*pow(phi,2));    
}
void DEF::system(double r,Array y,Array yprime,double beta)
{
    A = A_of_phi(y[4],beta);
    alpha = a_of_phi(y[4],beta);
    if(y[3]<eos.p_surface){eps=0.0;}
    else{eps=eos.e_at_p(y[3]);}
    yprime[0] = 4*M_PI*pow(A,4)*pow(r,2.0)*eps + 0.5*r*(r-2*y[0])*pow(y[2],2);
    yprime[1] = 2.0*(y[0] + 4.0*M_PI*pow(A,4)*pow(r,3.0)*y[3])/( r*(r-2.0*y[0])) + r*pow(y[2],2);
    yprime[2] =  4*M_PI*pow(A,4)*r*(alpha*(eps-3.0*y[3]) + r*(eps-y[3])*y[2] )/(r-2.0*y[0]) - 2*y[2]*(1 - y[0]/r)/(r - 2.0*y[0]);
    yprime[3] = -(eps + y[3])*((y[0] + 4.0*M_PI*pow(A,4)*pow(r,3.0)*y[3])/(r*(r-2.0*y[0])) + 0.5*r*pow(y[2],2) + alpha*y[2]);
    yprime[4] = y[2]; 

}

double DEF::dnudr_inf(double y0,double y2,double phi, double coupling)
{
    return 2.0*y0/(r_max*(r_max-2.0*y0))+r_max*pow(y2,2);
}
void DEF::computeMass(double& mass,double y0,double y2,double phi, double coupling)
{
    double dn_inf = dnudr_inf(y0, y2, phi,coupling);
    mass = dn_inf*pow(r_max,2)/(2*(1+r_max*dn_inf));
    // mass = dn_inf*pow(r_max,2)*exp(nu_inf)/2;
}
void DEF::computeRadius(double& radius, double r_surface, double phi_surface, double coupling)
{
    radius = r_surface*A_of_phi(phi_surface,coupling)*CGS::Length/pow(10,5);
}
void DEF::computeScalarCharge(double& scalarCharge, double dphi_last)
{
    scalarCharge = -dphi_last*pow(r_max,2);
}

R2::R2(std::string eos_name, double alpha,double foo):atol(1e-20),rtol(1e-06),reqmin(1e-15),
r_max(set_rmax(alpha)),r_max_surface(9.0),RDIV(5001),startingPointX(-1.0),startingPointY(0.1),
stepX(1.0),stepY(0.1), eos(eos_name){};

R2::~R2(){};

double R2::set_rmax(double alpha)
{
            if(alpha>=1.0 && alpha<20.0)
            {
                return 9.5;
            }
            else if(alpha>=pow(10,3))
            {
                return 500.0;
            }
            else
            {
                return alpha/2.0;
            }
}   
double R2::a_of_phi(double phi, double alpha)
{
    return -1/pow(3,0.5);
}
double R2::A_of_phi(double phi, double alpha)
{
    return exp(a_of_phi(phi,alpha)*phi);
}
double R2::V_of_phi(double phi, double alpha)
{
    return (1/(4*alpha))*pow(1 - exp(2*a_of_phi(phi,alpha)*phi),2);
}
double R2::dV_of_phi(double phi, double alpha)
{
    return (-a_of_phi(phi,alpha)/alpha)*(exp(a_of_phi(phi,alpha)*2*phi) - exp(a_of_phi(phi,alpha)*4*phi));
}
void R2::system(double r,Array y,Array yprime,double alpha)
{
    A = A_of_phi(y[4],alpha);
    V = V_of_phi(y[4],alpha);
    dV = dV_of_phi(y[4],alpha);
    if(y[3]<eos.p_surface)
        eps=0.0;
    else
        eps = eos.e_at_p(y[3]);
    yprime[0] = 4*M_PI*pow(A,4)*pow(r,2.0)*eps + 0.5*r*(r-2*y[0])*pow(y[2],2) + 0.25*V*pow(r,2);
    yprime[1] = 2.0*(y[0] + 4.0*M_PI*pow(A,4)*pow(r,3.0)*y[3] - 0.25*V*pow(r,3.0))/( r*(r-2.0*y[0])) + r*pow(y[2],2);
    yprime[2] = 4*M_PI*pow(A,4)*r*(R2::a_of_phi(y[4],alpha)*(eps-3.0*y[3]) + r*(eps-y[3])*y[2] )/(r-2.0*y[0]) - 2*y[2]*(1 - y[0]/r)/(r - 2.0*y[0]) + (0.5*V*y[2]*pow(r,2.0) + 0.25*dV*r)/(r - 2.0*y[0]);
    yprime[3] = -(eps + y[3])*((y[0] + 4.0*M_PI*pow(A,4)*pow(r,3.0)*y[3] - 0.25*V*pow(r,3.0))/(r*(r-2.0*y[0])) + 0.5*r*pow(y[2],2) + R2::a_of_phi(y[4],alpha)*y[2]);
    yprime[4] = y[2];
}
double R2::dnudr_inf(double y0,double y2,double phi, double coupling)
{
    double V = V_of_phi(phi,coupling);
    return 2.0*(y0 - 0.25*V*pow(r_max,3.0))/(r_max*(r_max-2.0*y0))+r_max*pow(y2,2);
}
void R2::computeMass(double& mass,double y0,double y2,double phi, double coupling)
{
    double dn_inf = dnudr_inf(y0, y2, phi,coupling);
    mass = dn_inf*pow(r_max,2)/(2*(1+r_max*dn_inf));
    // mass = dn_inf*pow(r_max,2)*exp(nu_inf)/2;
}
void R2::computeRadius(double& radius, double r_surface, double phi_surface, double coupling){
    radius = r_surface*A_of_phi(phi_surface,coupling)*CGS::Length/pow(10,5);
}
void R2::computeScalarCharge(double& scalarCharge, double dphi_last)
{
    scalarCharge = -dphi_last*pow(r_max,2);
}