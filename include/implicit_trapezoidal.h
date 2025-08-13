#ifndef IMPLICIT_TRAPEZOIDAL_H
#define IMPLICIT_TRAPEZOIDAL_H

#include "methods.h"

typedef struct {
    double y_n;
    double t_n;
    double t_next;
    double h;
} ImplicitTrapezoidalParams;

Result implicit_trapezoidal(double y0, double x0, double h, int n);

#endif