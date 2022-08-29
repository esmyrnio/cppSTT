%module pySTT
%{
	void get_MR(std::string theory,std::string eos_name, double central_density,double coupling,
    	double* mass, double* radius, double* scalarCharge, double* minimizationError, double* centralScalar);
%}

%include "std_string.i"
%include "typemaps.i"
%apply double *OUTPUT {double *mass, double *radius, double *scalarCharge, double* minimizationError,double* centralScalar};
%include <std_vector.i>


void get_MR(std::string theory,std::string eos_name, double central_density,double coupling,
double* mass, double* radius, double* scalarCharge, double* minimizationError, double* centralScalar);






