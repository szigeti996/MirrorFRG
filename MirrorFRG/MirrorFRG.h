// MirrorFRG.h : Include file for standard system include files,
// or project specific include files.

#pragma once

#include <iostream>
#include <cmath>
#include <array> 
#include <string>

const int n = 320;
const int n2 = 160;
const double pi = 3.14159265359;
const double mw = 783.;
const double EPS = 0.005;
const double size = 15.;

const double m0 = 185.; // MeV
const double muc = 922.7; // MeV
const double h1 = 17.59905821;
const double h2 = 9.953896925;
const double g = 9.45884097;

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

//template<typename T> // Bose-Einstein distribution at finite temperature
//T BoseEinstein(T omegaB){ return 1./(std::exp(omegaB/T0)-1);  }

//template<typename T> //Fermi-Dirac Temperature at finite temperature
//T FermiDirac(T omegaF, T mu){ return 1./(std::exp((omegaF-mu)/T0)-1.); }

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

double deriv(int j, int size, double step, std::array<double, n> f) { //Numerical derivation

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

    else if (j == 1) result = (f[2] - f[0]) / (2.0 * step);

    else if (j == size - 2) result = (f[size - 1] - f[size - 3]) / (2.0 * step);

    else if (j == 2) result = (f[0] - 8.0 * f[1] + 8.0 * f[3] - f[4]) / (12.0 * step);

    else if (j == size - 3) result = (f[size - 5] - 8.0 * f[size - 4] + 8.0 * f[size - 2] - f[size - 1]) / (12.0 * step);

    else result = (-f[j - 3] + 9.0 * f[j - 2] - 45.0 * f[j - 1] + 45.0 * f[j + 1] - 9.0 * f[j + 2] + f[j + 3]) / (60.0 * step);

    return result;
}

double deriv2(int j, int size, double step, std::array<double, n2> f) { //Numerical derivation

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

    else if (j == 1) result = (f[2] - f[0]) / (2.0 * step);

    else if (j == size - 2) result = (f[size - 1] - f[size - 3]) / (2.0 * step);

    else if (j == 2) result = (f[0] - 8.0 * f[1] + 8.0 * f[3] - f[4]) / (12.0 * step);

    else if (j == size - 3) result = (f[size - 5] - 8.0 * f[size - 4] + 8.0 * f[size - 2] - f[size - 1]) / (12.0 * step);

    else result = (-f[j - 3] + 9.0 * f[j - 2] - 45.0 * f[j - 1] + 45.0 * f[j + 1] - 9.0 * f[j + 2] + f[j + 3]) / (60.0 * step);

    return result;
}

std::array<double, n> func(double k, std::array<double, n> rho, std::array<double, n> Potential, std::array<double, n> mu) { //R.H.S of the Wetterich-flow equation
    std::array<double, n> UDer1;
    std::array<double, n> UDer2;
    std::array<double, n> Result;
    //std::array<double,n2> tempPot;
    //std::array<double,n2> tempDer1;
    //std::array<double,n2> tempDer2;

    for (int i = 0; i < n; i++)
        UDer1[i] = deriv(i, n, size, Potential);

    for (int i = 0; i < n; i++)
        UDer2[i] = deriv(i, n, size, UDer1);
    /*
    for(int i = 0; i < n2; i++)
        tempPot[i] = Potential[2*i];

    for(int i = 0; i < n2; i++){
        UDer1[2*i] = deriv2(i,n2,2.*size,tempPot);
        UDer1[2*i+1] = deriv2(i,n2,2.*size,tempPot);
        }
    */

    std::ofstream myfile2("pelda.dat");

    for (int i = 0; i < n; i++) {
        double Bosonterm = (1 + 2 * BoseEinsteinT0(Esigma(k, UDer1[i], UDer2[i], rho[i]))) / Esigma(k, UDer1[i], UDer2[i], rho[i]) + 3. * (1 + 2 * BoseEinsteinT0(Epi(k, UDer1[i]))) / Epi(k, UDer1[i]);
        double Nucleon1 = 4. / Enucleon(k, mass1(rho[i])) * (1 - FermiDiracT0(Enucleon(k, mass1(rho[i])), mu[i]) - FermiDiracT0(Enucleon(k, mass1(rho[i])), -mu[i]));
        double Nucleon2 = 4. / Enucleon(k, mass2(rho[i])) * (1 - FermiDiracT0(Enucleon(k, mass2(rho[i])), mu[i]) - FermiDiracT0(Enucleon(k, mass2(rho[i])), -mu[i]));
        myfile2 << rho[i] << "," << UDer1[i] << "," << UDer2[i] << "," << Esigma(k, UDer1[i], UDer2[i], rho[i]) << "," << Epi(k, UDer1[i]) << "," << std::pow(k, 4) / 12. / sqr(pi) * (Bosonterm - Nucleon1 - Nucleon2) << std::endl;
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

