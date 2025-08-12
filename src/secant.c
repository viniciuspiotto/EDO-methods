#include "secant.h"
#include <math.h>

double secant(function_with_params func, void* params, double x0, double x1, double epsilon, double delta) {
    double x2 = x1;
    int zero_division_error = 0;

    while (fabs(x2 - x1) > delta || fabs(func(x2, params)) > epsilon) {
        double f1 = func(x1, params);
        double f0 = func(x0, params);

        if (fabs(f1 - f0) == 0) {
            zero_division_error = 1;
            break;
        }

        x2 = x1 - f1 * (x1 - x0) / (f1 - f0);
        x0 = x1;
        x1 = x2;
    }

    return zero_division_error ? x1 : x2;
}