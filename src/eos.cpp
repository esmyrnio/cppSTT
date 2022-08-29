       
#include "../include/eos.hpp"
#include "../include/constants.hpp"
#include <cmath>
#include "../boost/math/interpolators/pchip.hpp"

using namespace SURFACE;
using namespace UNITS;
using namespace CGS;
using boost::math::interpolators::pchip;

EOS::EOS(std::string eos_file)
{
  std::string eos;
  if(EXECUTION::call_type=="main_call")
  {
    eos="EoS/";
  }
  else
  {
    eos="../EoS/";
  }
  eos+=eos_file;
  if((f_eos=fopen(eos.c_str(),"r")) == NULL ) {    
      std::cout<<"cannot open file "<<eos_file<<std::endl; 
      exit(0);
  }
  fscanf(f_eos,"%d\n",&n_tab);
  for(i=1;i<=n_tab;i++)
  {  
    fscanf(f_eos,"%lf %lf %lf %lf\n",&rho,&p,&h,&n0) ; 
    log_e_tab.push_back(log10(rho/Density));      /* multiply by C^2 to get */ 
    log_p_tab.push_back(log10(p/(Density*pow(C,2))));           /* energy density. */
    log_h_tab.push_back(log10(h/(C*C)));        
    log_n0_tab.push_back(log10(n0));                  /* STILL IN CGS ! */
  }
  p_surface=pow(10,log_p_tab[0]);
  e_surface=pow(10,log_e_tab[0]);
}

EOS::~EOS(){};

double EOS::e_at_p(double pp)
{
  if(pp<p_surface)
  {
    return 0;
  }
  else
  {
    auto spline = pchip<decltype(log_p_tab)>(log_p_tab,log_e_tab);
    return pow(10.0,spline(log10(pp)));
  }
}
double EOS::rho_at_p(double pp)
{
  if(pp<p_surface)
  {
    return 0;
  }
  else
  {
    auto spline = pchip<decltype(log_p_tab)>(log_p_tab,log_rho_tab);
    return pow(10.0,spline(log10(pp)));
  }
}
double EOS::p_at_e(double ee)
{
    auto spline = pchip<decltype(log_e_tab)>(log_e_tab,log_p_tab);
    return pow(10.0,spline(log10(ee)));
}
double EOS::p_at_h(double hh)
{ 
  if(hh<=enthalpy_min)
  {
    return 0;
  }
  else
  {
    auto spline = pchip<decltype(log_h_tab)>(log_h_tab,log_p_tab);
    return pow(10.0,spline(log10(hh)));
  }
}
double EOS::h_at_p(double pp)
{ 
    auto spline = pchip<decltype(log_p_tab)>(log_p_tab,log_h_tab);
    return pow(10.0,spline(log10(pp)));
}