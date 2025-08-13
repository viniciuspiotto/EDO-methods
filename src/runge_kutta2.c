#include "runge_kutta2.h"
#include "ode_function.h"
#include <stdlib.h>

Result runge_kutta2(double y0, double x0, double h, int n){
    int size = n + 1;
    Result res = {
        .size = size,
        .x_values = malloc(size * sizeof(double)),
        .y_values = malloc(size * sizeof(double))
    };

    res.x_values[0] = x0;
    res.y_values[0] = y0;

    for (int i = 1; i < res.size; i++){
        double x_prev = res.x_values[i-1];
        double y_prev = res.y_values[i-1];

        double k1 = ode_function(x_prev, y_prev);
        double k2 = ode_function(x_prev + h/2, y_prev + (h/2) * k1);

        res.x_values[i] = x_prev + h;
        res.y_values[i] = y_prev + h * k2;
    }

    return res;
}