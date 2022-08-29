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

using Vector = std::vector<double>;
using String = std::string;
using Array = double*;


template<class theoryType>
class STT{

    public:
        STT(String, double, double);
        ~STT();
    private:
        theoryType* theory;
        theoryType model;
        Vector linspace(double, double, int);
        Vector gradient(Vector, double);
        double functionToMinimize(double*);
        void findMin();
        void solve();
        typename TYPES<theoryType>::ODE_SYSTEM system;
    private:
        double eps,coupling, central_density, central_pressure, central_scalar, central_metric;
        double min_vals[2];
        double r,rout,dr;
        double m_c,nu_c,p_c;
        double m_1,nu_1,p_1;
        Vector r_vec;
        Vector y,yout,m,nu,phi,phi_trial,pr,dphi,dphi_trial,dnudr;
        int idx,idxlast,neq,istate;
        Vector r_s,m_s,p_s,e_s,rho_s,h_s,dmdr_s;
        double A, alpha;
        Vector startingPoints,steps;
        String eos_name;
    public:
        double mass,radius,scalarCharge,minimizationError,centralScalar;
        void computeModel();
        void printModelX();
        void printModelY();

};

template<class theoryType>
STT<theoryType>::STT(String eos_name, double central_density, double coupling):
        central_density(central_density*pow(10,15)/CGS::Density),coupling(coupling),
        model(eos_name,coupling,central_density),neq(5),
        central_pressure(model.eos.p_at_e(this->central_density)),
        theory(&model),
        eos_name(eos_name){}

template<class theoryType>
STT<theoryType>::~STT(){};

template<class theoryType>
Vector STT<theoryType>::linspace(double start, double end, int num)
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
template<class theoryType>
Vector STT<theoryType>::gradient(Vector input, double h){
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

template<class theoryType>
double STT<theoryType>::functionToMinimize(double* trial_vals){
  central_metric = trial_vals[0];
  central_scalar = trial_vals[1];
  // std::cout<<trial_vals[0]<<" "<<trial_vals[1]<<std::endl;
  r_vec = linspace(0.0,theory->r_max_surface, theory->RDIV);
  dr = r_vec[1]-r_vec[0];
  y = {0.0,central_metric,0.0,central_pressure,central_scalar};
  r = dr;
  rout = r + dr;
  LSODA<theoryType> lsoda;
  istate=1;

  idx=1;

  while(y[3]>1e-13 && r<r_vec[theory->RDIV-1]){
  lsoda.lsoda_update(theory, &theoryType::system, neq, y, yout, &r, rout, &istate,coupling,theory->rtol,theory->atol);
    rout += dr;
    y[0] = yout[1];
    y[1] = yout[2];
    y[2] = yout[3];
    y[3] = yout[4];
    y[4] = yout[5];
    idx += 1;
    if(istate==-3 || isnan(y[0]) || isnan(y[3])) break;
  }

  // std::cout<<" THIS IS r "<<r<<std::endl;

  idxlast=idx-1;
  lsoda.lsoda_update(theory, &theoryType::system , neq, y, yout, &r, theory->r_max, &istate,coupling,theory->rtol,theory->atol);
  rout += dr;
  y[0] = yout[1];
  y[1] = yout[2];
  y[2] = yout[3];
  y[3] = yout[4];
  y[4] = yout[5];
  idx += 1;
  // std::cout<<" THIS IS r2 "<<r<<std::endl;

  return pow(y[1],2)+pow(y[4],2);
}
template<class theoryType>
void STT<theoryType>::solve(){
    central_metric = min_vals[0];
    central_scalar = min_vals[1];
    r_vec = linspace(0.0,theory->r_max_surface, theory->RDIV);
    dr = r_vec[1]-r_vec[0];
    y = {0.0,central_metric,0.0,central_pressure,central_scalar};
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
    while(y[3]>1e-13 && r<r_vec[theory->RDIV-1]){
      m.push_back(y[0]);
      nu.push_back(y[1]);
      dphi.push_back(y[2]);
      pr.push_back(y[3]);
      phi.push_back(y[4]);
      lsoda.lsoda_update(theory, &theoryType::system, neq, y, yout, &r, rout, &istate,coupling,theory->rtol,theory->atol);
      rout += dr;
      y[0] = yout[1];
      y[1] = yout[2];
      y[2] = yout[3];
      y[3] = yout[4];
      y[4] = yout[5];
      idx += 1;
      if(istate==-3 || isnan(m[idx-2]) || isnan(pr[idx-2])) break;
  }
  idxlast=idx-1;
  lsoda.lsoda_update(theory, &theoryType::system , neq, y, yout, &r, theory->r_max, &istate,coupling,theory->rtol,theory->atol);
  rout += dr;
  y[0] = yout[1];
  y[1] = yout[2];
  y[2] = yout[3];
  y[3] = yout[4];
  y[4] = yout[5];
  idx += 1;
  m.push_back(y[0]);
  nu.push_back(y[1]);
  dphi.push_back(y[2]);
  pr.push_back(y[3]);
  phi.push_back(y[4]);
  
}
template<class theoryType>
void STT<theoryType>::findMin(){
    int number_of_variables = 2;
    double xmin[2],ynewlo;
    double start[2] = {model.startingPointX,model.startingPointY};
    double steps[2] = {model.stepX,model.stepY};
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
template<class theoryType>
void STT<theoryType>::computeModel(){
    findMin();
    solve();
    theory->computeMass(mass,m.back(),dphi.back(),phi.back(),coupling);
    theory->computeRadius(radius,r_vec[idxlast],phi[idxlast],coupling);
    theory->computeScalarCharge(scalarCharge,dphi.back());
}

template<class theoryType>
void STT<theoryType>::printModelX(){
  printf("\n");
  printf("  %s                  \n\n",eos_name.c_str());
  printf("  %2.3e     CENTRAL DENSITY         (10^15 gr/cm^3)\n",central_density*CGS::Density/pow(10,15));
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
void STT<theoryType>::printModelY(){
  printf("\n");
  printf("  %s                  \n\n",eos_name.c_str());
  printf("  %2.3e     CENTRAL DENSITY         (10^15 gr/cm^3)\n",central_density*CGS::Density/pow(10,15));
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