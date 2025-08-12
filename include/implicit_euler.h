#ifndef IMPLICIT_EULER_H
#define IMPLICIT_EULER_H

#include "methods.h"

typedef struct {
    double y_n;
    double t_next;
    double h;
} EulerImplicitParams;

Result implicit_euler(double y0, double x0, double h, int n);

#endif