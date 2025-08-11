#ifndef METHODS_H
#define METHODS_H

typedef struct result {
    double * x_values;
    double * y_values;
    int size;
} Result;

extern double g_y_n;
extern double g_t_next;
extern double g_h;
extern double g_y_n;
extern double g_t_n;

double f(double x0, double y0);
Result euler (double y0, double x0, double h, double n);
Result adamsBashford2 (double y0, double x0, double h, double n);
Result adamsMoulton2 (double y0, double y1, double x0, double h, double n);
Result BDF2 (double y0, double y1, double x0, double h, double n);
Result eulerImplicito (double y0, double x0, double h, double n);
Result trapezioImplicito (double y0, double x0, double h, double n);

#endif