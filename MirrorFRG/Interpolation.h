#pragma once
#include <iostream>
#include <cmath>
#include <gsl/gsl_spline.h>

struct my_f_params { gsl_interp_accel* facc; gsl_spline* fspline; };

double function1(double x, void* p) {
    struct my_f_params* params = (struct my_f_params*)p;
    gsl_interp_accel* facc = (params->facc);
    gsl_spline* fspline = (params->fspline);

    double f = gsl_spline_eval(fspline, x, facc);
    return f;
}