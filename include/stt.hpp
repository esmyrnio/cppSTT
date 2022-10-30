#ifndef STT_HPP
#define STT_HPP
#include <math.h>
#include <vector>
#include <stdlib.h>
#include "LSODA.h"
#include "asa047.hpp"
#include "eos.hpp"
#include "constants.hpp"
#include "theory.hpp"
#include "types.hpp"
#include <stdexcept>

using String = std::string;

/* This is the main class, which is templated w.r.t gravitational theory type (DEF or R2).The modified TOV system is solved here. */

template<class theoryType>
class STT{

    public:
        STT(String, double, double); //constructor taking the path to EoS file, central energy density and theory's coupling
        ~STT();

    private:
        Vector linspace(double, double, int) const; //linear space fector
        Vector gradient(Vector const&, double) const; //find gradient of vector
        double functionToMinimize(const Array); //returns boundary conditions at infinity 
        void findMin(); //this is the double shooting method routine for boundary conditions at infinity 
        void solve(); //solves the ODE system

    private:
        theoryType* theory;
        theoryType model;
        double coupling, centralDensity, centralPressure, centralMetric;
        double r,rout,dr; //LSODA integration points and grid's step size
        Vector r_vec; 
        Vector y,yout,m,nu,phi,pr,dphi; //solution vectors
        int idx,idxlast,neq,istate; //  index and LSODA parameters
        double min_vals[2]; // values obtained from the Nelder-Mead optimization routine
        Vector steps; // steps for the Nelder-Mead optimization routine
        String eos_name; // path to EoS file
        typename TYPES<theoryType>::ODE_SYSTEM system; // system of ODEs
        
    public:
        double mass, radius, scalarCharge, minimizationError, centralScalar;
        void computeModel();
        void printModelX() const;
        void printModelY() const;
};


template<class theoryType>
STT<theoryType>::STT(String eos_name, double centralDensity, double coupling):
        centralDensity(centralDensity*pow(10,15)/CGS::Density),coupling(coupling),
        model(eos_name,coupling,centralDensity),neq(5),
        centralPressure(model.eos.p_at_e(this->centralDensity)),
        theory(&model),
        eos_name(eos_name){}

template<class theoryType>
STT<theoryType>::~STT(){};


// returns a linear space grid starting from "star" to "end" with "num" grid points 
template<class theoryType>
Vector STT<theoryType>::linspace(double start, double end, int num) const
{
  Vector linspaced;
  if (num == 0) { return linspaced; }
  if (num == 1) 
    {
      linspaced.push_back(start);
      return linspaced;
    }
  double delta = (end - start) / (num - 1);
  for(int i=0; i < num-1; ++i)
    {
      linspaced.push_back(start + delta * i);
    }
  linspaced.push_back(end); 
  return linspaced;
}

// returns the gradient of a vector
template<class theoryType>
Vector STT<theoryType>::gradient(Vector const& input, double h) const
{
    if (input.size() <= 1) return input;
    Vector res;
    for(int j=0; j<input.size(); j++) {
        int j_left = j - 1;
        int j_right = j + 1;
        if (j_left < 0) {
            j_left = 0;
            j_right = 1;
        }
        if (j_right >= input.size()){
            j_right = input.size() - 1; 
            j_left = j_right - 1;
        }
        double dist_grad = (input[j_right] - input[j_left]) / (2.0*h);
        res.push_back(dist_grad);
    }
    return res;
}

// returns the metric and scalar value at infinity to be optimized to zero
template<class theoryType>
double STT<theoryType>::functionToMinimize(const Array trial_vals)
{  
  centralMetric = trial_vals[0];
  centralScalar = trial_vals[1];

  r_vec = linspace(0.0,theory->r_max_surface, theory->RDIV);
  dr = r_vec[1]-r_vec[0];
  y = {0.0,centralMetric,0.0,centralPressure,centralScalar};
  r = dr;
  rout = r + dr;
  LSODA<theoryType> lsoda;
  istate=1;
  idx=1;
  
  // integration stops when pressure drops below given threshold
  while(y[3]>1e-13 && r<r_vec[theory->RDIV-1])
  {
    lsoda.lsoda_update(theory, &theoryType::system, neq, y, yout, &r, rout, &istate,coupling,theory->rtol,theory->atol);
    rout += dr;
    y = {yout[1],yout[2],yout[3],yout[4],yout[5]};
    idx += 1;
    if(istate==-3 || isnan(y[0]) || isnan(y[3])) break;
  }
  
  idxlast=idx-1;
  lsoda.lsoda_update(theory, &theoryType::system , neq, y, yout, &r, theory->r_max, &istate,coupling,theory->rtol,theory->atol);
  rout += dr;
  y = {yout[1],yout[2],yout[3],yout[4],yout[5]};
  idx += 1;
  return pow(y[1],2)+pow(y[4],2);
}

// Nelder-Mead routine, filling the min_vals array
template<class theoryType>
void STT<theoryType>::findMin()
{
    int number_of_variables = 2;
    double xmin[2],ynewlo;
    double start[2] = {model.metricStartingPoint,model.scalarStartingPoint};
    double steps[2] = {model.metricStep,model.scalarStep};
    int icount,konvge,kcount,numres,ifault;
    konvge = 10;
    kcount = 500;
    nelmin(this, &STT<theoryType>::functionToMinimize, number_of_variables, start, xmin, &ynewlo,
    model.reqmin, steps, konvge, kcount, &icount, &numres, &ifault);
    min_vals[0] = xmin[0];
    min_vals[1] = xmin[1];
    minimizationError = ynewlo;
    centralScalar = xmin[1];
}

// main routine
template<class theoryType>
void STT<theoryType>::solve()
{
  centralMetric = min_vals[0];
  centralScalar = min_vals[1];
  r_vec = linspace(0.0,theory->r_max_surface, theory->RDIV);
  dr = r_vec[1]-r_vec[0];
  y = {0.0,centralMetric,0.0,centralPressure,centralScalar};
  m.push_back(y[0]);
  nu.push_back(y[1]);
  dphi.push_back(y[2]);
  pr.push_back(y[3]);
  phi.push_back(y[4]);
  istate=1;
  r = dr;
  rout = r + dr;
  idx=1;
  LSODA<theoryType> lsoda;

  while(y[3]>1e-13 && r<r_vec[theory->RDIV-1])
  {
    m.push_back(y[0]);
    nu.push_back(y[1]);
    dphi.push_back(y[2]);
    pr.push_back(y[3]);
    phi.push_back(y[4]);
    lsoda.lsoda_update(theory, &theoryType::system, neq, y, yout, &r, rout, &istate,coupling,theory->rtol,theory->atol);
    rout += dr;
    y = {yout[1],yout[2],yout[3],yout[4],yout[5]};
    idx += 1;
    if(istate==-3 || isnan(m[idx-2]) || isnan(pr[idx-2])) break;
  }

  idxlast=idx-1;
  lsoda.lsoda_update(theory, &theoryType::system , neq, y, yout, &r, theory->r_max, &istate,coupling,theory->rtol,theory->atol);
  rout += dr;
  y = {yout[1],yout[2],yout[3],yout[4],yout[5]};
  idx += 1;
  m.push_back(y[0]);
  nu.push_back(y[1]);
  dphi.push_back(y[2]);
  pr.push_back(y[3]);
  phi.push_back(y[4]);
}

// computes model
template<class theoryType>
void STT<theoryType>::computeModel()
{
    findMin();
    solve();
    theory->computeMass(mass,m.back(),dphi.back(),phi.back(),coupling);
    theory->computeRadius(radius,r_vec[idxlast],phi[idxlast],coupling);
    theory->computeScalarCharge(scalarCharge,dphi.back());
}

template<class theoryType>
void STT<theoryType>::printModelX() const
{
  printf("\n");
  printf("  %s                  \n\n",eos_name.c_str());
  printf("  %2.3e     CENTRAL DENSITY         (10^15 gr/cm^3)\n",centralDensity*CGS::Density/pow(10,15));
  printf("\n");
  printf("  %2.2e     COUPLING                (km^2)\n",coupling);
  printf("\n");
  printf("  %2.3e     BOUNDARY MINIMIZATION             \n",std::abs(minimizationError));
  printf("  %2.2e     SCALAR AT CENTER             \n",centralScalar);
  printf("\n");
  printf("  %2.3e     MASS                    (M_sun)\n",mass);
  printf("  %2.3e     RADIUS                  (km)\n",radius);
  printf("  %2.2e     SCALAR CHARGE                  \n",scalarCharge);
  printf("\n");
}

template<class theoryType>
void STT<theoryType>::printModelY() const
{
  printf("\n");
  printf("  %s                  \n\n",eos_name.c_str());
  printf("  %2.3e     CENTRAL DENSITY         (10^15 gr/cm^3)\n",centralDensity*CGS::Density/pow(10,15));
  printf("\n");
  printf("  %2.3e     COUPLING                (km^2)\n",coupling);
  printf("\n");
  printf("  %2.3e     BOUNDARY MINIMIZATION             \n",std::abs(minimizationError));
  printf("  %2.3e     SCALAR AT CENTER             \n",centralScalar);
  printf("\n");
  printf("  %2.3e     MASS                    (M_sun)\n",mass);
  printf("  %2.3e     RADIUS                  (km)\n",radius);
  printf("  %2.3e     SCALAR CHARGE                  \n",scalarCharge);
  printf("\n");
}
#endif