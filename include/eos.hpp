#ifndef EOS_HPP
#define EOS_HPP

#include <iostream>
#include <vector>
#include <string>
#include <fstream>

class EOS{
    public:
        EOS(std::string);
        ~EOS();

        double p_at_e(double ee);
        double e_at_p(double pp);
        double p_at_h(double hh);
        double h_at_p(double pp);
        double rho_at_p(double pp);

    public:
        double e_surface,p_surface;
    private:
        int i,n_tab;
        double p,rho,h,n0;
        FILE* f_eos;
        std::vector<double> log_e_tab,log_p_tab,log_h_tab,log_n0_tab,log_rho_tab;

};

#endif