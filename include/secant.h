#ifndef SECANT_H
#define SECANT_H

typedef double (*function_with_params)(double, void*);

double secant(function_with_params func, void* params, double x0, double x1, double epsilon, double delta);

#endif