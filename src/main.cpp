#include "../include/stt.hpp"
#include <iostream>
#include <stdlib.h>
#include <string.h>

// function to be exposed to python using SWIG. Computes model's mass, radius and scalar charge.
void get_MR(std::string theory,std::string eos_name, double central_density,double coupling,
            double* mass, double* radius, double* scalarCharge,
            double* minimizationError, double* centralScalar)
{
    if(theory=="DEF")
    {
        STT<DEF> model(eos_name,central_density,coupling);
        model.computeModel();
        *mass = model.mass;
        *radius = model.radius;
        *scalarCharge = model.scalarCharge;
        *minimizationError = model.minimizationError;
        *centralScalar = model.centralScalar;

    }
    else if(theory=="R2")
    {
        STT<R2> model(eos_name,central_density,coupling);
        model.computeModel();
        *mass = model.mass;
        *radius = model.radius;
        *scalarCharge = model.scalarCharge;
        *minimizationError = model.minimizationError;
        *centralScalar = model.centralScalar;
    }
}

int main(int argc, char *argv[]){
    char eos_name[50];
    char theory[50];
    double central_density, coupling;
    EXECUTION::CALL_TYPE = "main_call";
    for(int i=1;i<argc;i++) 
      if(argv[i][0]=='-'){
        switch(argv[i][1]){
          case 'm':
            sscanf(argv[i+1],"%s",theory);
            break;
          case 'f':
            sscanf(argv[i+1],"%s",eos_name);
            break;
          case 'e':
            sscanf(argv[i+1],"%lf",&central_density);
            break;
          case 'c':
            sscanf(argv[i+1],"%lf",&coupling);
            break;
          }          
      }
    if(!strcmp(theory,"DEF"))
    {
        STT<DEF> model(eos_name,central_density,coupling);
        model.computeModel();
        model.printModelX();
    }
    else if(!strcmp(theory,"R2"))
    {
        STT<R2> model(eos_name,central_density,coupling);
        model.computeModel();
        model.printModelY();
    }
    else
    {
        std::cout<<"theory not available."<<std::endl;
    }
    return 0;
}