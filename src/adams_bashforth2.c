#include "adams_bashforth2.h"
#include "ode_function.h"
#include <stdlib.h>

Result adams_bashforth2(double y0, double x0, double h, int n) {
    int size = n + 1;
    Result res = {
        .size = size,
        .x_values = malloc(size * sizeof(double)),
        .y_values = malloc(size * sizeof(double))
    };

    res.x_values[0] = x0;
    res.y_values[0] = y0;

    res.x_values[1] = x0 + h;
    res.y_values[1] = y0 + h * ode_function(x0, y0);

    for (int i = 1; i < n; i++) {
        res.x_values[i + 1] = res.x_values[i] + h;

        double f_current = ode_function(res.x_values[i], res.y_values[i]);
        double f_prev = ode_function(res.x_values[i - 1], res.y_values[i - 1]);

        res.y_values[i + 1] = res.y_values[i] + h * ((3.0 / 2.0) * f_current - (1.0 / 2.0) * f_prev);
    }

    return res;
}