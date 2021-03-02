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

const int n = 400;
const double pi = 3.14159265359;
const double mw = 783.;
const double EPS = 0.005;
const double size = 12.;

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

template<typename T> // Bose-Einstein distribution at zero temperature
T BoseEinsteinT0(T omegaB) {
    if (std::abs(omegaB) < 1.0e-5)
        return 1.;
    else
        return 0.;
}

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

std::array<double, n> func(double k, std::array<double, n> rho, std::array<double, n> Potential, std::array<double, n> mu) { //R.H.S of the Wetterich-flow equation
    
    std::array<double, n> Result;

    double* r; r = (double*)malloc((size_t)sizeof(double) * n);
    memset(r, 0, (size_t)sizeof(double) * n);
    double* Pot; Pot = (double*)malloc((size_t)sizeof(double) * n);
    memset(Pot, 0, (size_t)sizeof(double) * n);

    double* UDer1; UDer1 = (double*)malloc((size_t)sizeof(double) * n);
    memset(UDer1, 0, (size_t)sizeof(double) * n);
    double* UDer2; UDer2 = (double*)malloc((size_t)sizeof(double) * n);
    memset(UDer2, 0, (size_t)sizeof(double) * n);

    for (int i = 0; i < n; i++) {
        r[i] = rho[i];
        Pot[i] = Potential[i];
    }
    
    gsl_interp_accel* facc = gsl_interp_accel_alloc();
    gsl_spline* fspline = gsl_spline_alloc(gsl_interp_cspline, n);
    gsl_spline_init(fspline, r, Pot, n);
    struct my_f_params params = { facc, fspline };
    gsl_function F;
    F.function = &function1;
    F.params = &params;

    double result1, abserr1;
    gsl_deriv_forward(&F, r[0], 1., &result1, &abserr1);
    UDer1[0] = result1;
    
    for (int i = 1; i < n; i++) {
        if (i == n-1) {
            double result, abserr;
            gsl_deriv_backward(&F, r[i], 0.01, &result, &abserr);
            UDer1[i] = result;
        }
        else {
            double result, abserr;
            gsl_deriv_central(&F, r[i], 0.01, &result, &abserr);
            UDer1[i] = result;
        }
    }
    
    gsl_spline_init(fspline, r, UDer1, n);
    struct my_f_params params1 = { facc, fspline };
    gsl_function G;
    G.function = &function1;
    G.params = &params1;
    
    gsl_deriv_forward(&G, r[0], 1., &result1, &abserr1);
    UDer2[0] = result1;

    for (int i = 1; i < n; i++) {
        if (i == n - 1) {
            double result, abserr;
            gsl_deriv_backward(&G, r[i], 0.01, &result, &abserr);
            UDer2[i] = result;
        }
        else {
            double result, abserr;
            gsl_deriv_central(&G, r[i], 0.01, &result, &abserr);
            UDer2[i] = result;
        }
    }
   
    const int N = 101;

    double* r2; r2 = (double*)malloc((size_t)sizeof(double) * N);
    memset(r2, 0, (size_t)sizeof(double) * N);
    double* Es; Es = (double*)malloc((size_t)sizeof(double) * N);
    memset(Es, 0, (size_t)sizeof(double) * N);

    for (int i = 0; i < N-1; i++) {
        int j = i * 4;
        r2[i] = r[j];
        Es[i] = Esigma(k, UDer1[j], UDer2[j], r[j]);
    }

    r2[N - 1] = r[n - 1];
    Es[N - 1] = Esigma(k, UDer1[n - 1], UDer2[n - 1], r[n - 1]);

    gsl_interp_accel* facc1 = gsl_interp_accel_alloc();
    gsl_spline* fspline1 = gsl_spline_alloc(gsl_interp_cspline, N);
    gsl_spline_init(fspline1, r2, Es, N);
    struct my_f_params params2 = { facc1, fspline1 };
    gsl_function H;
    H.function = &function1;
    H.params = &params2;

    std::ofstream myfile2("pelda.dat");
    for (int i = 0; i < n; i++) {
        double Sigma = GSL_FN_EVAL(&H, r[i]);
        double Bosonterm = (1 + 2 * BoseEinsteinT0(Sigma)) / Sigma + 3. * (1 + 2 * BoseEinsteinT0(Epi(k, UDer1[i]))) / Epi(k, UDer1[i]);
        double Nucleon1 = 4. / Enucleon(k, mass1(rho[i])) * (1 - FermiDiracT0(Enucleon(k, mass1(rho[i])), mu[i]) - FermiDiracT0(Enucleon(k, mass1(rho[i])), -mu[i]));
        double Nucleon2 = 4. / Enucleon(k, mass2(rho[i])) * (1 - FermiDiracT0(Enucleon(k, mass2(rho[i])), mu[i]) - FermiDiracT0(Enucleon(k, mass2(rho[i])), -mu[i]));
        myfile2 << rho[i] << " " << UDer1[i] << " " << UDer2[i] << " " << Sigma << " " << Epi(k, UDer1[i]) << " " << std::pow(k, 4) / 12. / sqr(pi) * (Bosonterm - Nucleon1 - Nucleon2) << std::endl;
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

