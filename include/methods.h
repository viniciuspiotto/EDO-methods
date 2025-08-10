#ifndef METHODS_H
#define METHODS_H

#define PI 3.14159265358979323846

typedef struct result {
    double * x_values;
    double * y_values;
    int size;
} Result;

Result euler (double y0, double x0, double h, double n);
Result adamsBashford2 (double y0, double x0, double h, double n);
Result adamMulton2 (double y0, double x0, double h, double n);
Result BDF2 (double y0, double x0, double h, double n);
Result eulerImplicito (double y0, double x0, double h, double n);
Result trapezioImplicito (double y0, double x0, double h, double n);

#endif