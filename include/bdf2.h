#ifndef BDF2_H
#define BDF2_H

#include "methods.h"

typedef struct {
    double y_n_prev;
    double y_n;
    double t_next;
    double h;
} BDF2Params;

Result bdf2(double y0, double y1, double x0, double h, int n);

#endif