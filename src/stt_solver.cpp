#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>
#include <iomanip>

#include "../include/LSODA.h"

#include "../include/asa047.h"

#define G 6.6732e-8
#define C 2.9979e+10
#define MSUN 1.989e+33
#define PI 3.1415926535
#define MB 1.66e-24
#define N 45001  // grid size

// #define __USE_MINGW_ANSI_STDIO  (__MINGW_FEATURES__ & __MINGW_ANSI_STDIO__)


char eos_name[80];

double Length = G*MSUN/pow(C,2),
Time = Length/C,
Density = MSUN/pow(Length,3),
central_density,
log_e_tab[201],               /* rho points in tabulated EOS */
log_p_tab[201],               /* p points in tabulated EOS */
log_h_tab[201],              /* h points in EOS file */
log_n0_tab[201],         /* number density in EOS file */  
log_rho_tab[201],
scalar_charge,
Mass;                                // gravitational mass


double coupling;

int n_tab,                           /* number of tabulated EOS points */
    n_nearest=1,                     /* nearest grid point, used in interp. */ 
    print_option;                    /* select print out */ 

void load_eos( char eos_file[])
{
 std::size_t i;                    /* counter */

 double p,                 /* pressure */
        rho,               /* density */
        h,                 /* enthalpy */
        n0;                /* number density */    
      
 FILE *f_eos;              /* pointer to eos_file */
 

    /* OPEN FILE TO READ */
    char eos[100];
    strcat(eos,"EoS/");  /* add path */
    strcat(eos,eos_file);
	
    if((f_eos=fopen(eos,"r")) == NULL ) {    
       printf("cannot open file:  %s\n",eos_file); 
       exit(0);
    }

 
    /* READ NUMBER OF TABULATED POINTS */

    fscanf(f_eos,"%d\n",&n_tab);


    /* READ EOS, H, N0 AND MAKE THEM DIMENSIONLESS (EXCEPT N0) */
 
    for(i=1;i<=n_tab;i++) {  
      /*  input line for _four_column equation of state file */
       fscanf(f_eos,"%lf %lf %lf %lf\n",&rho,&p,&h,&n0) ; 

/*  input line for _five_ column equation of state file         
    fscanf(f_eos,"%lf %lf %lf %lf %*lf\n",&rho,&p,&h,&n0) ; */

       log_e_tab[i]=log10(rho/Density);       /* multiply by C^2 to get */ 
       log_p_tab[i]=log10(p/(Density*pow(C,2)));             /* energy density. */
       log_h_tab[i]=log10(h/(C*C));        
       log_n0_tab[i]=log10(n0);
       log_rho_tab[i] = log10(MB*n0*pow(Length,3)/MSUN);             /* STILL IN CGS ! */
    }
}

void hunt(double xx[], int n, double x, int *jlo)
{ 
  int jm,jhi,inc,ascnd;

  ascnd=(xx[n] > xx[1]);
  if (*jlo <= 0 || *jlo > n) {
    *jlo=0;
    jhi=n+1;
  } else {
    inc=1;
    if (x >= xx[*jlo] == ascnd) {
      if (*jlo == n) return;
      jhi=(*jlo)+1;
      while (x >= xx[jhi] == ascnd) {
        *jlo=jhi;
        inc += inc;
        jhi=(*jlo)+inc;
        if (jhi > n) {
          jhi=n+1;
          break;
        }
      }
    } else {
      if (*jlo == 1) {
        *jlo=0;
        return;
      }
      jhi=(*jlo);
      *jlo -= 1;
      while (x < xx[*jlo] == ascnd) {
        jhi=(*jlo);
        inc += inc;
        *jlo=jhi-inc;
        if (*jlo < 1) {
          *jlo=0;
          break;
        }
      }
    }
  }
  while (jhi-(*jlo) != 1) {
    jm=(jhi+(*jlo)) >> 1;
    if (x > xx[jm] == ascnd)
      *jlo=jm;
    else
      jhi=jm;
  }
}

/*******************************************************************/
double interp(double xp[], double yp[], int np ,double xb)
{ 
 int k,        /* index of 1st point */
     m=4;      /* degree of interpolation */ 
 
 double y;     /* intermediate value */

 hunt(xp,np,xb,&n_nearest);

 k=std::min(std::max(n_nearest-(m-1)/2,1),np+1-m);

 if( xb==xp[k] ||  xb==xp[k+1] || xb==xp[k+2] || xb==xp[k+3]) xb += 1e-12;

 y= (xb-xp[k+1])*(xb-xp[k+2])*(xb-xp[k+3])*yp[k]/
        ((xp[k]-xp[k+1])*(xp[k]-xp[k+2])*(xp[k]-xp[k+3]))
 
    +(xb-xp[k])*(xb-xp[k+2])*(xb-xp[k+3])*yp[k+1]/
       ((xp[k+1]-xp[k])*(xp[k+1]-xp[k+2])*(xp[k+1]-xp[k+3]))
 
    +(xb-xp[k])*(xb-xp[k+1])*(xb-xp[k+3])*yp[k+2]/
       ((xp[k+2]-xp[k])*(xp[k+2]-xp[k+1])*(xp[k+2]-xp[k+3]))
 
    +(xb-xp[k])*(xb-xp[k+1])*(xb-xp[k+2])*yp[k+3]/
       ((xp[k+3]-xp[k])*(xp[k+3]-xp[k+1])*(xp[k+3]-xp[k+2]));

 return (y);
}
/*******************************************************************/
/*******************************************************************/

double e_at_p(double pp)
{
   return pow(10.0,interp(log_p_tab,log_e_tab,n_tab,log10(pp)));
}
/*******************************************************************/
double e_at_rho(double rhorho)
{
   return pow(10.0,interp(log_rho_tab,log_e_tab,n_tab,log10(rhorho))); 
}
/*******************************************************************/
double p_at_e(double ee)
{
  return pow(10.0,interp(log_e_tab,log_p_tab,n_tab,log10(ee)));
}
/*******************************************************************/
double p_at_rho(double rhorho)
{
  return pow(10.0,interp(log_rho_tab,log_p_tab,n_tab,log10(rhorho)));
}
/*******************************************************************/
double p_at_h(double hh)
{ 
  return pow(10.0,interp(log_h_tab,log_p_tab,n_tab,log10(hh)));
}
/*******************************************************************/
double h_at_p(double pp)
{ 
  return pow(10.0,interp(log_p_tab,log_h_tab,n_tab,log10(pp)));
}
/*******************************************************************/
double n0_at_e(double ee)
{ 
  return pow(10.0,interp(log_e_tab,log_n0_tab,n_tab,log10(ee)));
}
/*******************************************************************/
/*******************************************************************/
/*******************************************************************/

static double A_of_phi(double phi,double beta){
	return exp(0.5*beta*pow(phi,2));
}
static double alpha_of_phi(double phi, double beta){
        return beta*phi;
}

/*******************************************************************/
/*******************************************************************/
/*******************************************************************/


/*******************************************************************/
/*******************************************************************/

template<typename T>
std::vector<double> linspace(T start_in, T end_in, int num_in)
{

  std::vector<double> linspaced;

  double start = static_cast<double>(start_in);
  double end = static_cast<double>(end_in);
  double num = static_cast<double>(num_in);

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

/*******************************************************************/
/*******************************************************************/


/**************************************************************************************************************************************/

static void system_inside(double r, double* y, double* yprime,double beta)
{
    
    (void)r;

    double A = A_of_phi(y[4],beta);

    double alpha = alpha_of_phi(y[4],beta);

    double eps = e_at_p(y[3]);


    yprime[0] = 4*PI*pow(A,4)*pow(r,2.0)*eps + 0.5*r*(r-2*y[0])*pow(y[2],2);

   
    yprime[1] = 2.0*(y[0] + 4.0*PI*pow(A,4)*pow(r,3.0)*y[3])/( r*(r-2.0*y[0])) + r*pow(y[2],2);

 	  yprime[2] =  4*PI*pow(A,4)*r*(alpha*(eps-3.0*y[3]) + r*(eps-y[3])*y[2] )/(r-2.0*y[0]) - 2*y[2]*(1 - y[0]/r)/(r - 2.0*y[0]);

    yprime[3] = -(eps + y[3])*((y[0] + 4.0*PI*pow(A,4)*pow(r,3.0)*y[3])/(r*(r-2.0*y[0])) + 0.5*r*pow(y[2],2) + alpha*y[2]);

    yprime[4] = y[2];
            
}

/**************************************************************************************************************************************/

static void system_outside(double r, double* y, double* yprime,double beta)
{
    
    (void)r;

    yprime[0] = 0.5*r*(r-2*y[0])*pow(y[2],2);

    yprime[1] = 2.0*(y[0])/(r*(r-2.0*y[0])) + r*pow(y[2],2);
    
    yprime[2] = -2.0*y[2]*(1 - y[0]/r)/(r - 2.0*y[0]);

    yprime[3] = 0;

    yprime[4] = y[2];
            
}

/**************************************************************************************************************************************/
/**************************************************************************************************************************************/

std::vector<double> gradient(std::vector<double> input, double h){
    if (input.size() <= 1) return input;
    std::vector<double> res;
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
        // gradient value at position j
        double dist_grad = (input[j_right] - input[j_left]) / (2.0*h);
        res.push_back(dist_grad);
    }
    return res;
}

/**************************************************************************************************************************************/
/**************************************************************************************************************************************/

/*double phi_c_peak(double veta)
{
  return -1.0174 - 4.0407/veta;
}

double phi_c_width(double veta)
{
  return 0.0015*pow(0.2684,veta);
}

double phi_c_guess(double e_c, double veta)
{
  double phi;
  phi = phi_c_width(veta)*pow(e_c-1.5,2) + phi_c_peak(veta);

  if(phi>=0){
    return 0;
  }
  else{
    return phi;
  }
}*/

/**************************************************************************************************************************************/
/**************************************************************************************************************************************/

double function_to_minimize(double args[2])
{
  double met_phi_c_obj,
  phi_c_obj,
  eps_c_obj,
  p_c_obj,
  r_max_obj;

  met_phi_c_obj = args[0];
  phi_c_obj = args[1];


  eps_c_obj = central_density*pow(10,15)/Density;
  p_c_obj = p_at_e(eps_c_obj); 

	
	r_max_obj = 128*4.0 *pow(3.0/(4.0*PI*eps_c_obj),1.0/3.0);


  int neq = 5;
  int istate_obj = 1;
  double t_obj, tout_obj;

  std::vector<double> r_vec_obj = linspace(0.0,r_max_obj,N);
  double dr_obj = r_vec_obj[1]-r_vec_obj[0];

  std::vector<double> y_obj = {0.0,met_phi_c_obj,0.0,p_at_e(eps_c_obj),phi_c_obj};


	t_obj = dr_obj;
	tout_obj = t_obj+dr_obj;

  LSODA lsoda;
  

  std::vector<double> yout_obj,
    yyout_obj,
    res_la_obj,
    res_ph_obj,
    res_om_obj,
    res_pr_obj,
    res_phi_obj;


  res_la_obj.push_back(y_obj[0]);
  res_ph_obj.push_back(y_obj[1]);
  res_om_obj.push_back(y_obj[2]);
  res_pr_obj.push_back(y_obj[3]);
  res_phi_obj.push_back(y_obj[4]);


  int idx_obj=2;


  while(y_obj[3]>=1e-10 && t_obj<r_vec_obj[N-1]){

  	lsoda.lsoda_update(system_inside, neq, y_obj, yout_obj, &t_obj, tout_obj, &istate_obj, coupling,1e-06,1e-20);
      res_la_obj.push_back(yout_obj[1]);
      res_ph_obj.push_back(yout_obj[2]);
      res_om_obj.push_back(yout_obj[3]);
      res_pr_obj.push_back(yout_obj[4]);
      res_phi_obj.push_back(yout_obj[5]);


      if(istate_obj==-3 || isnan(res_ph_obj[idx_obj-2]) || isnan(res_phi_obj[idx_obj-2])){
        break;
      }
      
      tout_obj += dr_obj;

      y_obj[0] = yout_obj[1];
      y_obj[1] = yout_obj[2];
      y_obj[2] = yout_obj[3];
      y_obj[3] = yout_obj[4];
      y_obj[4] = yout_obj[5];

      idx_obj += 1;
  }


  int idxlast_obj = idx_obj-1;


  while(t_obj<r_vec_obj[N-1]){

	lsoda.lsoda_update(system_outside, 5, y_obj, yout_obj, &t_obj, tout_obj, &istate_obj, coupling,1e-06,1e-20);
      res_la_obj.push_back(yout_obj[1]);
      res_ph_obj.push_back(yout_obj[2]);
      res_om_obj.push_back(yout_obj[3]);
      res_pr_obj.push_back(yout_obj[4]);
      res_phi_obj.push_back(yout_obj[5]);

      if(istate_obj==-3 || isnan(res_ph_obj[idx_obj-2]) || isnan(res_phi_obj[idx_obj-2])){
        break;
      }
      
      tout_obj += dr_obj;

      y_obj[0] = yout_obj[1];
      y_obj[1] = yout_obj[2];
      y_obj[2] = yout_obj[3];
      y_obj[3] = yout_obj[4];
      y_obj[4] = yout_obj[5];

      idx_obj += 1;    	

  }

return pow(pow(res_ph_obj[idx_obj-3],2)+pow(res_phi_obj[idx_obj-3],2),0.5);
}

/**************************************************************************************************************************************/
/**************************************************************************************************************************************/

double *find_min(double e_c, double b)
{

  static double min_vals[2];

  int number_of_variables = 2;

  double starting_points[2];

  double xmin[2];
  double ynewlo;
  double step[2];

  int icount;
  int konvge,kcount;
  int numres;
  int ifault;

  double reqmin = 1.0E-12;

  starting_points[0] = -1.0;
  starting_points[1] = -0.5;
  //starting_points[1] = phi_c_guess(e_c,b);

  step[0] = 1.0;
  step[1] = -0.1;

  konvge = 10;
  kcount = 500;
  
  nelmin (function_to_minimize, number_of_variables, starting_points, xmin, &ynewlo, reqmin, step,
    konvge, kcount, &icount, &numres, &ifault);


  if (abs(ynewlo)>1e-06){


    double starting_points2[2];

    starting_points2[0] = -1.0;
    starting_points2[1] = -0.5;

    double xmin2[2];
    double ynewlo2;
    double step2[2];

    double reqmin2 = 1.0E-12;
  
    int icount2;
    int konvge2,kcount2;
    int numres2;
    int ifault2;

    step2[0] = 1.0;
    step2[1] = 0.1;

    konvge2 = 10;
    kcount2 = 500;

    nelmin (function_to_minimize, 2, starting_points2, xmin2, &ynewlo2, reqmin2, step2,
    konvge2, kcount2, &icount2, &numres2, &ifault2);


    std::cout << "\n";
    std::cout << "\n";

    std::cout << "  Estimate of minimizing value X*:\n";
    std::cout << "\n";
    for (std::size_t i = 0; i < 2; i++ )
    {
      if (i==0)
      {
        std::cout << "   "<<"metric_phi_central = " << setw(7) << xmin2[i] << "\n";
      }
      else{
        std::cout << "   "<<"scalar_central = " << setw(7) << xmin2[i] << "\n";
      }
      
    }

    std::cout << "\n";
    std::cout << "\n";
    std::cout << "  F(X*) = " << ynewlo2 << "\n";

    std::cout << "\n";
    std::cout << "  Number of Nelder-Mead iterations = " << icount2 << "\n";
    std::cout << "  Number of Nelder-Mead restarts =   " << numres2 << "\n";
    std::cout << "\n";
    std::cout << "\n";

    min_vals[0] = xmin2[0];
    min_vals[1] = xmin2[1];
    printf("-------------------------------------------------------------------------------\n");  
  
  }

  else{

    std::cout << "\n";
    std::cout << "\n";

    std::cout << "  Estimate of minimizing value X*:\n";
    std::cout << "\n";
    for (std::size_t i = 0; i < 2; i++ )
    {
      if (i==0)
      {
        std::cout << "   "<<"metric_phi_central = " << setw(7) << xmin[i] << "\n";
      }
      else{
        std::cout << "   "<<"scalar_central = " << setw(7) << xmin[i] << "\n";
      }
      
    }

    std::cout << "\n";
    std::cout << "\n";
    std::cout << "  F(X*) = " << ynewlo << "\n";

    std::cout << "\n";
    std::cout << "  Number of Nelder-Mead iterations = " << icount << "\n";
    std::cout << "  Number of Nelder-Mead restarts =   " << numres << "\n";
    std::cout << "\n";
    std::cout << "\n";

    min_vals[0] = xmin[0];
    min_vals[1] = xmin[1];
    printf("-------------------------------------------------------------------------------\n");  

  }

  return min_vals;
}

/**************************************************************************************************************************************/
/**************************************************************************************************************************************/
/**************************************************************************************************************************************/
/**************************************************************************************************************************************/

int main(int argc, char const *argv[])
{

  for(int i=1;i<argc;i++) 
      if(argv[i][0]=='-'){
        switch(argv[i][1]){
          
          case 'f':
              sscanf(argv[i+1],"%s",eos_name);
              break;
          case 'c':
              sscanf(argv[i+1],"%lf",&coupling);
              break;
          case 'e':
              sscanf(argv[i+1],"%lf",&central_density);
              break;
          case 'p':
              sscanf(argv[i+1],"%i",&print_option);
              break;
          }
      }

  load_eos(eos_name);

  double *min_val;
	
  min_val = find_min(central_density,coupling);

  double met_phi_c,
        phi_c,
        eps_c,
        p_c,
        r_max;
  
  met_phi_c = min_val[0];
  phi_c = min_val[1];

  eps_c = central_density*pow(10,15)/Density;
  p_c = p_at_e(eps_c); 
  
  r_max = 128*4.0 *pow(3.0/(4.0*PI*eps_c),1.0/3.0);

  int neq = 5;
  int istate = 1;
  double t, tout;

  std::vector<double> r_vec = linspace(0.0,r_max,N);
  std::vector<double> y = {0.0,met_phi_c,0.0,p_at_e(eps_c),phi_c};

  double dr = r_vec[1]-r_vec[0];
  double met_lamda_c,omega_c,pr_c;

  t = dr;
  tout = t+dr;

  LSODA lsoda;

  std::vector<double> yout,
                 yyout,
                 res_la,
                 res_ph,
                 res_om,
                 res_pr,
                 res_phi;


  res_la.push_back(y[0]);
  res_ph.push_back(y[1]);
  res_om.push_back(y[2]);
  res_pr.push_back(y[3]);
  res_phi.push_back(y[4]);

  int idx=2;

  while(y[3]>=1e-10 && t<r_vec[N-1]){

    lsoda.lsoda_update(system_inside, 5, y, yout, &t, tout, &istate, coupling,1e-06,1e-20);
      res_la.push_back(yout[1]);
      res_ph.push_back(yout[2]);
      res_om.push_back(yout[3]);
      res_pr.push_back(yout[4]);
      res_phi.push_back(yout[5]);
      
      if(istate==-3 || isnan(res_ph[idx-2]) || isnan(res_phi[idx-2])){
      break;
    }
      tout += dr;

      y[0] = yout[1];
      y[1] = yout[2];
      y[2] = yout[3];
      y[3] = yout[4];
      y[4] = yout[5];

      idx += 1;

  }

  int idxlast = idx-1;

  while(t<r_vec[N-1]){

  lsoda.lsoda_update(system_outside, 5, y, yout, &t, tout, &istate, coupling,1e-06,1e-20);
      res_la.push_back(yout[1]);
      res_ph.push_back(yout[2]);
      res_om.push_back(yout[3]);
      res_pr.push_back(yout[4]);
      res_phi.push_back(yout[5]);

      if(istate==-3 || isnan(res_ph[idx-2]) || isnan(res_phi[idx-2])){
      break;
    }
      
      tout += dr;

      y[0] = yout[1];
      y[1] = yout[2];
      y[2] = yout[3];
      y[3] = yout[4];
      y[4] = yout[5];

      idx += 1;     
  }

  std::vector<double> dphdr;
  std::vector<double> dphidr;

  dphdr = gradient(res_ph,dr);
  dphidr = gradient(res_phi,dr);

  Mass = dphdr[N-3]*pow(r_vec[N-3],2)*exp(res_ph[N-3])/2;
        
  //Mass = dphdr[N-3]*pow(r_vec[N-3],2)/(2*(1+r_vec[N-3]*dphdr[N-3]));
  
  scalar_charge = -dphidr[N-3]*pow(r_vec[N-3],2);

  switch(print_option){

    case 0:

    printf("\n");
    printf("  %s                  \n",eos_name);    
    printf("  %i         N             \n",N);
    printf("  %6.5e  e_c           (10^15 gr/cm^3)\n",central_density);
    printf("\n");
    printf(" %2.2e     beta             \n",coupling);
    printf("\n");
    printf("  %6.5e  M             (M_sun)\n",Mass);
    printf("  %6.5e  R             (km)\n",r_vec[idxlast]*A_of_phi(res_phi[idxlast],coupling)*Length/pow(10,5));
    printf(" %6.5e  scalar charge          \n",scalar_charge);
    break;

    case 1:

    printf("\n");
    printf("  %s                  \n",eos_name);
    printf("  %i         N             \n",N);
    printf("  %6.5e  e_c           (10^15 gr/cm^3)\n",central_density);
    printf("\n");
    printf(" %2.2e     beta             \n",coupling);
    printf("\n");
    printf("  %6.5e  M             (M_sun)\n",Mass);
    printf("  %6.5e  R             (km)\n",r_vec[idxlast]*A_of_phi(res_phi[idxlast],coupling)*Length/pow(10,5));
    printf(" %6.5e  scalar charge          \n",scalar_charge);
    printf("\n");


    printf("-------------------------------------------------------------------------------\n");   
    printf("r        mu         nu        phi        epsilon        pressure\n");
    printf("-------------------------------------------------------------------------------\n");   


    for(std::size_t w =0;w<=N;w++){

        printf("%5.4e %5.4e %5.4e %5.4e %5.4e %5.4e \n",
          r_vec[w], res_la[w], res_ph[w], res_phi[w], e_at_p(res_pr[w]), res_pr[w]);
        }

    break;
  }

  return 0;
}
