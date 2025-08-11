#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include "methods.h"

double g_y_n_prev;
double g_y_n;
double g_t_next;
double g_h;

double f (double x, double y){
    return 10 * (12 * cos(2 * M_PI * x) - y);
}

double gEulerImplicito (double y) {
    return y - g_y_n - g_h * f(g_t_next, y);
}

Result euler (double y0, double x0, double h, double n){
    Result res;
    res.size = n + 1;

    res.x_values = malloc(res.size * sizeof(double));
    res.y_values = malloc(res.size * sizeof(double));

    res.x_values[0] = x0;
    res.y_values[0] = y0;

    for (int i = 1; i <= n; i++){
        res.x_values[i] = res.x_values[i-1] + h;
        res.y_values[i] = res.y_values[i-1] + h * f(res.x_values[i-1], res.y_values[i-1]);
    }

    return res;
}

double secante(double (*f)(double),double x0, double x1, double epsilon, double delta) {
    double x2 = x1;  
    int zeroDivisionError = 0;

    while (fabs(x2 - x1) > delta || fabs(f(x2)) > epsilon) {
        if (fabs(f(x1) - f(x0)) == 0) {
            zeroDivisionError = 1;
            break;
        } else {
            x2 = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0));
        }
        x0 = x1;
        x1 = x2;
    }

    if (zeroDivisionError) {
        return x1;
    }

    return x2;
}

Result eulerImplicito (double y0, double x0, double h, double n){

    Result res;
    res.size = n + 1;

    res.x_values = malloc(res.size * sizeof(double));
    res.y_values = malloc(res.size * sizeof(double));

    res.x_values[0] = x0;
    res.y_values[0] = y0;

    for (int i = 1; i <= n; i++){
        res.x_values[i] = res.x_values[i-1] + h;

        g_y_n = res.y_values[i-1];
        g_t_next = res.x_values[i];
        g_h = h;

        res.y_values[i] = secante(gEulerImplicito, g_y_n, g_y_n + 0.1, 1e-8, 1e-8);
    }

    return res;
}

double gBDF2(double y) {
    return 3*y - 2*g_h*f(g_t_next, y) - (4*g_y_n-g_y_n_prev);
}

Result BDF2(double y0, double y1, double x0, double h, double n) {
    Result res;
    res.size = n + 1;

    res.x_values = malloc(res.size * sizeof(double));
    res.y_values = malloc(res.size * sizeof(double));

    res.x_values[0] = x0;
    res.y_values[0] = y0;
    res.x_values[1] = x0 + h;
    res.y_values[1] = y1;

    for (int i = 2; i <= n; i++) {
        res.x_values[i] = res.x_values[i-1] + h;

        g_y_n_prev = res.y_values[i-2];
        g_y_n = res.y_values[i-1];
        g_t_next = res.x_values[i];
        g_h = h;

        res.y_values[i] = secante(gBDF2, g_y_n, g_y_n + 0.1, 1e-8, 1e-8);
    }

    return res;
}