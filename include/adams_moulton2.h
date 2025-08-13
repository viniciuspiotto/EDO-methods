#ifndef ADAMS_MOULTON2_H
#define ADAMS_MOULTON2_H

#include "methods.h"

typedef struct {
    double y_n;
    double t_next;
    double h;
    double f_k;
    double f_k_1;
} AdamsMoulton2Params;

Result adams_moulton2(double y0, double y1, double x0, double h, int n);

#endif