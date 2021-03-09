// MirrorFRG.cpp : Defines the entry point for the application.
#include <iostream> 
#include <fstream>
#include <vector>
#include <cmath>
#include <array> 
#include <string>

#include "MirrorFRG.h"

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_spline.h>

int main()
{
    double k = 360.; //MeV
    double Lambda = 361.; //MeV
    double dk = 1000 * 0.000001; // MeV
    double w0start = -90.; //Mev

    std::array<double, n> rho;
    std::array<double, n> U0;
    std::array<double, n> mu;

    std::vector<double> tmp1, tmp2, tmp3;
    std::vector<double> tmp4, tmp5, tmp6;

    std::ifstream infile("test.dat");
    double a, b, c;
    while (infile >> a >> b >> c) {
        tmp1.push_back(a);
        tmp2.push_back(b);
        tmp3.push_back(c);
    }
    
    for (int i = 0; i < n - 1; i++) {
        int j = i * 50;
        rho[i] = tmp1[j];
        mu[i] = tmp3[j];
        U0[i] = tmp2[j];
    }

    rho[n - 1] = tmp1[2199];
    mu[n - 1] = tmp3[2199];
    U0[n - 1] = tmp2[2199];

    PotentialT0 U = { rho,U0,mu };
    
    tmp1.clear();
    tmp2.clear();
    tmp3.clear();

    for (int i = 0; i < 1; i++){
        U.Omega0Rho(U, k, w0start);
        U.Euler(U, k, dk);
        k = k + dk;
        std::cout << k << std::endl;
    }

    std::ofstream myfile("Result.dat");
    for (int i = 0; i < n; i++) {
        myfile << U.rho[i] << " " << U.Uk[i] << " " << U.mu[i] << std::endl;
    }

    myfile.close();
    return 0;
}
