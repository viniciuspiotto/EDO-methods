#ifndef IMPLICIT_MIDPOINT_H
#define IMPLICIT_MIDPOINTT_H

#include "methods.h"

typedef struct {
    double y_n;       
    double t_n;       
    double h;         
} ImplicitMidpointParams;

Result implicit_midpoint(double y0, double x0, double h, int n);

#endif