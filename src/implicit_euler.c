#include "implicit_euler.h"
#include "ode_function.h"
#include "secant.h"
#include <stdlib.h>

static double implicit_euler_equation(double y, void* params) {
    EulerImplicitParams* p = (EulerImplicitParams*) params;
    return y - p->y_n - p->h * ode_function(p->t_next, y);
}

Result implicit_euler(double y0, double x0, double h, int n) {
    int size = n + 1;
    Result res = {
        .size = size,
        .x_values = malloc(size * sizeof(double)),
        .y_values = malloc(size * sizeof(double))
    };

    res.x_values[0] = x0;
    res.y_values[0] = y0;

    EulerImplicitParams params = { .h = h };
    for (int i = 1; i <= n; i++) {
        res.x_values[i] = res.x_values[i - 1] + h;
        params.y_n = res.y_values[i - 1];
        params.t_next = res.x_values[i];
        res.y_values[i] = secant(implicit_euler_equation, &params, params.y_n, params.y_n + 0.1, 1e-8, 1e-8);
    }

    return res;
}