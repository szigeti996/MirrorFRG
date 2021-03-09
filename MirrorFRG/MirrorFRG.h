// MirrorFRG.h : Include file for standard system include files,
// or project specific include files.

#include <iostream>
#include <cmath>
#include <array> 
#include <string>

#include "Interpolation.h"

#include <gsl/gsl_spline.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>

#define GSL_FN_EVAL(F,x) (*((F)->function))(x,(F)->params)

const int n = 45;

const double pi = 3.14159265359;
const double mw = 783.;
const double EPS = 0.005;

const double m0 = 185.; // MeV
const double muc = 922.7; // MeV
const double h1 = 17.59905821;
const double h2 = 9.953896925;
const double g = 9.45884097;

//template<typename T> // Bose-Einstein distribution at finite temperature
//T BoseEinstein(T omegaB){ return 1./(std::exp(omegaB/T0)-1);  }

//template<typename T> //Fermi-Dirac Temperature at finite temperature
//T FermiDirac(T omegaF, T mu){ return 1./(std::exp((omegaF-mu)/T0)-1.); }

template <typename T> //square
T sqr(T x) { return x * x; }

template <typename T> //postive mass (mirror)
T mass1(T rho) { return 1. / std::sqrt(2) * (std::sqrt(rho * sqr(h1 + h2) + sqr(m0)) + std::sqrt(rho) * (h1 - h2)); }

template <typename T> //negative mass (nucleon)
T mass2(T rho) { return 1. / std::sqrt(2) * (std::sqrt(rho * sqr(h1 + h2) + sqr(m0)) - std::sqrt(rho) * (h1 - h2)); }

template <typename T>  //Energy of the sigma meson
T Esigma(T k, T Uderiv, T Uderiv2, T rho) { return std::sqrt(k * k + 2 * Uderiv2 * rho + Uderiv); }

template <typename T> //Energy of the Pi meson
T Epi(T k, T Uderiv) { return std::sqrt(k * k + Uderiv); }

template<typename T> //Energy of the fermions
T Enucleon(T k, T mass) { return std::sqrt(k * k + mass * mass); }

template<typename T> //Fermi-Dirac Temperature at zero temperature
T FermiDiracT0(T omegaF, T mu) {
    if (std::abs(omegaF - mu) < 1.0e-5)
        return 1.;
    else
        return 0.;
}

template<typename T> //Dirac-delta function
T DiracDelta(T x, T x0) {
    if (std::abs(x - x0) < 1.0e-7)
        return true;
    else
        return false;
}

template<typename T> //Numerical evaluation of the integral contains Dirac-Delta
T w0func(T k0, T mu, T m) { return sqr(sqr(k0)) / 12. / sqr(pi) / std::sqrt(sqr(k0) + sqr(m)) / std::abs(std::sqrt(1 - sqr(m / mu))); }

double w0equation(double rho, double k, double w0) { //w0(rho,,k) equation at zero temperature
    double mueff = muc + g * w0;
    double m1 = mass1(rho);
    double m2 = mass2(rho);
    double k01 = sqr(mueff) - sqr(m1);
    double k02 = sqr(mueff) - sqr(m2);
    //std::cout << k01 << " " << k02 << std::endl;
    if (k01 > 0 && k02 > 0) {
        bool eq1 = DiracDelta(k, std::sqrt(k01));
        bool eq2 = DiracDelta(k, std::sqrt(k02));
        if (eq1 == true && eq2 == true)
            return sqr(mw) * w0 + g * w0func(std::sqrt(k01), mueff, m1) + g * w0func(std::sqrt(k02), mueff, m2);
        if (eq1 == true && eq2 == false)
            return sqr(mw) * w0 + g * w0func(std::sqrt(k01), mueff, m1);
        if (eq1 == false && eq2 == true)
            return sqr(mw) * w0 + g * w0func(std::sqrt(k02), mueff, m2);
        else
            return sqr(mw) * w0;
    }
    if (k01 > 0 && k02 < 0) {
        bool eq1 = DiracDelta(k, std::sqrt(k01));
        if (eq1 == true)
            return sqr(mw) * w0 + g * w0func(std::sqrt(k01), mueff, m1);
        else
            return sqr(mw) * w0;
    }
    if (k02 > 0 && k01 < 0) {
        bool eq2 = DiracDelta(k, std::sqrt(k02));
        if (eq2 == true)
            return sqr(mw) * w0 + g * w0func(std::sqrt(k02), mueff, m2);
        else
            return sqr(mw) * w0;
    }
    else
        return sqr(mw) * w0;
}

double Root_Finder(double rho, double k, double a, double b) {
    double c = a;
    while ((b - a) >= EPS) {
        // Find middle point
        c = (a + b) / 2;
        // Check if middle point is root
        if (w0equation(rho, k, c) == 0.0)
            break;
        // Decide the side to repeat the steps
        else if (w0equation(rho, k, c) * w0equation(rho, k, a) < 0)
            b = c;
        else
            a = c;
    }
    return c;
}

double deriv(int j, int size, double step, std::array<double, n> const & f) { //Numerical derivation

    double result;

    if (j == 0) {
        double A = (f[2] - f[0]) / (2.0 * step);
        double B = (f[0] - 8.0 * f[1] + 8.0 * f[3] - f[4]) / (12.0 * step);
        result = 2 * A - B;
    }

    else if (j == size - 1) {
        double A = (f[size - 1] - f[size - 3]) / (2.0 * step);
        double B = (f[size - 5] - 8.0 * f[size - 4] + 8.0 * f[size - 2] - f[size - 1]) / (12.0 * step);
        result = 2 * A - B;
    }

    else if(j == 1) 
        result = (f[2] - f[0]) / (2.0 * step);
    else if(j == 2)
        result = (f[0] - 8.0 * f[1] + 8.0 * f[3] - f[4]) / (12.0 * step);

    else if (j == size - 2)
        result = (f[size - 1] - f[size - 3]) / (2.0 * step);

    else if(j == size - 3)
        result = (f[size - 2] - f[size - 4]) / (2.0 * step);

    else if(j == size - 4)
        result = (f[size - 3] - f[size - 5]) / (2.0 * step);

    else if (j == size - 5)
        result = (f[size - 4] - f[size - 6]) / (2.0 * step);

    else
        result = (-f[j - 3] + 9.0 * f[j - 2] - 45.0 * f[j - 1] + 45.0 * f[j + 1] - 9.0 * f[j + 2] + f[j + 3]) / (60.0 * step);

    return result;
}

std::array<double, n> func(double k, std::array<double, n> const& rho, std::array<double, n> const& Potential, std::array<double, n> const& mu) { //R.H.S of the Wetterich-flow equation

    std::ofstream myfile2("pelda.dat");
    std::array<double, n> Result;
    std::array<double, n> Uder1;
    std::array<double, n> Uder2;

    const double size = 50.;

    for (int i = 0; i < n; i++)
        Uder1[i] = deriv(i, n, size, Potential);

    for (int i = 0; i < n; i++)
        Uder2[i] = deriv(i, n, size, Uder1);

    for (int i = 0; i < n; i++) {
        double Sigma = Esigma(k, Uder1[i], Uder2[i], rho[i]);
        double pion = Epi(k, Uder1[i]);
        double Bosonterm = 1. / Sigma + 3. / pion;
        double Nucleon1 = 4. / Enucleon(k, mass1(rho[i])) * (1 - FermiDiracT0(mu[i], Enucleon(k, mass1(rho[i]))) - FermiDiracT0(-mu[i], Enucleon(k, mass1(rho[i]))));
        double Nucleon2 = 4. / Enucleon(k, mass2(rho[i])) * (1 - FermiDiracT0(mu[i], Enucleon(k, mass2(rho[i]))) - FermiDiracT0(-mu[i], Enucleon(k, mass2(rho[i]))));
        myfile2 << rho[i] << " " << Potential[i] << " " << Uder1[i] << " " << Uder2[i] << " " << Sigma << " " << pion << " " << std::pow(k, 4) / 12. / sqr(pi) * (Bosonterm - Nucleon1 - Nucleon2) << std::endl;
        Result[i] = std::pow(k, 4) / 12. / sqr(pi) * (Bosonterm - Nucleon1 - Nucleon2);
    }
    myfile2.close();

    return Result;
}

struct PotentialT0 { //Potential function at every k;
    std::array<double, n> rho, Uk, mu;

    PotentialT0& Omega0Rho(PotentialT0 const& v, double k, double w0start) {
        for (int i = 0; i < n; i++) {
            double w0root = Root_Finder(v.rho[i], k, w0start, 0.);
            if (w0root > -0.005)
                w0root = 0.;
            mu[i] = v.mu[i] + g * w0root;
        }
        return *this;
    }

    PotentialT0& Euler(PotentialT0 const& v, double k, double dk) {

        std::array<double, n> tmp;
        tmp = func(k, v.rho, v.Uk, v.mu);
        for (int i = 0; i < n; i++)
            Uk[i] = v.Uk[i] + dk * tmp[i];
        return *this;
    }

};


/*   
    double* r; r = (double*)malloc((size_t)sizeof(double) * N);
    memset(r, 0, (size_t)sizeof(double) * N);
    double* der1; der1 = (double*)malloc((size_t)sizeof(double) * N);
    memset(der1, 0, (size_t)sizeof(double) * N);
    double* der2; der2 = (double*)malloc((size_t)sizeof(double) * N);
    memset(der2, 0, (size_t)sizeof(double) * N);

int size = 50;

    for (int i = 0; i < N - 1; i++) {
        int j = i * size;
        potent[i] = Potential[j];
        r[i] = rho[j];
    }

    potent[N - 1] = Potential[n - 1];
    r[N - 1] = rho[n - 1];

    for (int i = 0; i < N; i++)
        Uder1[i] = deriv(i, N, (double)size, potent);

    for (int i = 0; i < N; i++)
        Uder2[i] = deriv(i, N, (double)size, Uder1);

    for (int i = 0; i < N; i++){
        der1[i] = Uder1[i];
        der2[i] = Uder2[i];
    }
 
    for (int i = 0; i < N; i++)
        myfile2 << r[i] << " " << Uder1[i] << " " << Uder2[i] << " " << Esigma(k, Uder1[i], Uder2[i], rho[i]) <<std::endl;

    gsl_interp_accel* facc = gsl_interp_accel_alloc();
    gsl_spline* fspline = gsl_spline_alloc(gsl_interp_cspline_periodic, N);
    gsl_spline_init(fspline, r, der1, N);
    struct my_f_params params = { facc, fspline };
    gsl_function F;
    F.function = &function1;
    F.params = &params;

    for (int i = 0; i < n; i++)
        d1[i] = GSL_FN_EVAL(&F, rho[i]);

    gsl_spline_init(fspline, r, der2, N);
    struct my_f_params params2 = { facc, fspline };
    gsl_function G;
    G.function = &function1;
    G.params = &params2;*/


/*
 const int N = 250;

    double* r; r = (double*)malloc((size_t)sizeof(double) * N);
    memset(r, 0, (size_t)sizeof(double) * N);
    double* Pot; Pot = (double*)malloc((size_t)sizeof(double) * N);
    memset(Pot, 0, (size_t)sizeof(double) * N);
    double* UDer1; UDer1 = (double*)malloc((size_t)sizeof(double) * N);
    memset(UDer1, 0, (size_t)sizeof(double) * N);
    double* UDer2; UDer2 = (double*)malloc((size_t)sizeof(double) * N);
    memset(UDer2, 0, (size_t)sizeof(double) * N);

    for (int i = 0; i < N - 1; i++) {
        int j = i * 10;
        r[i] = rho[j];
        Pot[i] = Potential[j];
    }

    r[N - 1] = rho[n - 1];
    Pot[N - 1] = Potential[n - 1];

    gsl_interp_accel* facc = gsl_interp_accel_alloc();
    gsl_spline* fspline = gsl_spline_alloc(gsl_interp_cspline, N);

    gsl_spline_init(fspline, r, Pot, N);
    struct my_f_params params = { facc, fspline };
    gsl_function F;
    F.function = &function1;
    F.params = &params;

    double result1, abserr1;
    gsl_deriv_forward(&F, r[0], 1., &result1, &abserr1);
    UDer1[0] = result1;

    for (int i = 1; i < N; i++) {
        if (i == N-1) {
            double result, abserr;
            gsl_deriv_backward(&F, r[i], 1, &result, &abserr);
            UDer1[i] = result;
        }
        else {
            double result, abserr;
            gsl_deriv_central(&F, r[i], 0.01, &result, &abserr);
            UDer1[i] = result;
        }
    }

    gsl_spline_init(fspline, r, UDer1, N);
    struct my_f_params params1 = { facc, fspline };
    gsl_function G;
    G.function = &function1;
    G.params = &params1;

    for (int i = 0; i < n; i++)
        der1[i] = GSL_FN_EVAL(&G, rho[i]);

    gsl_deriv_forward(&G, r[0], 1., &result1, &abserr1);
    UDer2[0] = result1;

    for (int i = 1; i < N; i++) {
        if (i == N - 1) {
            double result, abserr;
            gsl_deriv_backward(&G, r[i], 1., &result, &abserr);
            UDer2[i] = result;
        }
        else {
            double result, abserr;
            gsl_deriv_central(&G, r[i], 1., &result, &abserr);
            UDer2[i] = result;
        }
    }

    gsl_spline_init(fspline, r, UDer2, N);
    struct my_f_params params2 = { facc, fspline };
    gsl_function H;
    H.function = &function1;
    H.params = &params2;
*/

/*std::array<double, n> func(double k,std::array<double, n> const & rho, std::array<double, n> const & Potential, std::array<double, n> const & mu) { //R.H.S of the Wetterich-flow equation
    
    std::ofstream myfile2("pelda.dat");
    std::array<double, n> Result;
    std::array<double, n> Uder1;
    std::array<double, n> Uder2;

    const double size = 50.; 

    for (int i = 0; i < n; i++)
        Uder1[i] = deriv(i, n, size, Potential);

    for (int i = 0; i < n; i++)
        Uder2[i] = deriv(i, n, size, Uder1);
    
    for (int i = 0; i < n; i++) {
        double Sigma = Esigma(k, Uder1[i], Uder2[i], rho[i]);
        double pion = Epi(k, Uder1[i]);
        double Bosonterm = 1./Sigma + 3. / pion;
        double Nucleon1 = 4. / Enucleon(k, mass1(rho[i])) * (1 - FermiDiracT0(mu[i],Enucleon(k, mass1(rho[i]))) - FermiDiracT0(-mu[i],Enucleon(k, mass1(rho[i]))));
        double Nucleon2 = 4. / Enucleon(k, mass2(rho[i])) * (1 - FermiDiracT0(mu[i],Enucleon(k, mass2(rho[i]))) - FermiDiracT0(-mu[i],Enucleon(k, mass2(rho[i]))));
        myfile2 << rho[i] << " " << Uder1[i] << " " << Uder2[i] << " " << Sigma << " " << pion << " " << std::pow(k, 4) / 12. / sqr(pi) * (Bosonterm - Nucleon1 - Nucleon2) << std::endl;
        Result[i] = std::pow(k, 4) / 12. / sqr(pi) * (Bosonterm - Nucleon1 - Nucleon2);
    }
    myfile2.close();
    
    return Result;
}
*/

/*
std::array<double, n> func(double k,std::array<double, n> const & rho, std::array<double, n> const & Potential, std::array<double, n> const & mu) { //R.H.S of the Wetterich-flow equation

    std::ofstream myfile2("pelda.dat");
    std::ofstream myfile3("pelda2.dat");
    std::array<double, n> Result;
    std::array<double, n> Es, Ep;
    const double step = 50.;

    std::array<double, N> potent;
    std::array<double, N> Uder1, Uder2;

    double* r; r = (double*)malloc((size_t)sizeof(double) * N);
    memset(r, 0, (size_t)sizeof(double) * N);
    double* Esig; Esig = (double*)malloc((size_t)sizeof(double) * N);
    memset(Esig, 0, (size_t)sizeof(double) * N);
    double* Epion; Epion = (double*)malloc((size_t)sizeof(double) * N);
    memset(Epion, 0, (size_t)sizeof(double) * N);

    for (int i = 0; i < N - 1; i++) {
        potent[i] = Potential[i * 50];
        r[i] = rho[i * 50];
    }

    r[N - 1] = rho[n - 1];
    potent[N - 1] = Potential[n - 1];

    for (int i = 0; i < N; i++)
        Uder1[i] = deriv(i, N, step, potent);

    for (int i = 0; i < N; i++)
        Uder2[i] = deriv(i, N, step, Uder1);

    for (int i = 0; i < N; i++) {
        Esig[i] = Esigma(k, Uder1[i], Uder2[i], r[i]);
        myfile3 << r[i] << " " << Esig[i] << std::endl;
    }
    for (int i = 0; i < N; i++)
        Epion[i] = Epi(k, Uder1[i]);

    gsl_interp_accel* facc = gsl_interp_accel_alloc();
    gsl_spline* fspline = gsl_spline_alloc(gsl_interp_cspline_periodic, N);
    gsl_spline_init(fspline, r, Esig, N);
    struct my_f_params params = { facc, fspline };
    gsl_function F;
    F.function = &function1;
    F.params = &params;

    for (int i = 0; i < n; i++)
        Es[i] = GSL_FN_EVAL(&F, rho[i]);

    gsl_spline_init(fspline, r, Epion, N);
    struct my_f_params params1 = { facc, fspline };
    gsl_function G;
    G.function = &function1;
    G.params = &params1;

    for (int i = 0; i < n; i++)
        Ep[i] = GSL_FN_EVAL(&G, rho[i]);

    for (int i = 0; i < n; i++) {
        double Bosonterm = 1./Es[i] + 3. / Ep[i];
        double Nucleon1 = 4. / Enucleon(k, mass1(rho[i])) * (1 - FermiDiracT0(mu[i],Enucleon(k, mass1(rho[i]))) - FermiDiracT0(-mu[i],Enucleon(k, mass1(rho[i]))));
        double Nucleon2 = 4. / Enucleon(k, mass2(rho[i])) * (1 - FermiDiracT0(mu[i],Enucleon(k, mass2(rho[i]))) - FermiDiracT0(-mu[i],Enucleon(k, mass2(rho[i]))));
        myfile2 << rho[i] << " " << Es[i] << " " << Ep[i] << " " << std::pow(k, 4) / 12. / sqr(pi) * (Bosonterm - Nucleon1 - Nucleon2) << std::endl;
        Result[i] = std::pow(k, 4) / 12. / sqr(pi) * (Bosonterm - Nucleon1 - Nucleon2);
    }
    myfile2.close();

    return Result;
}
*/