#include "bdf2.h"
#include "ode_function.h"
#include "secant.h"
#include <stdlib.h>

static double bdf2_equation(double y, void* params) {
    BDF2Params* p = (BDF2Params*) params;
    return 3 * y - 2 * p->h * ode_function(p->t_next, y) - (4 * p->y_n - p->y_n_prev);
}

Result bdf2(double y0, double y1, double x0, double h, int n) {
    int size = n + 1;
    Result res = {
        .size = size,
        .x_values = malloc(size * sizeof(double)),
        .y_values = malloc(size * sizeof(double))
    };

    res.x_values[0] = x0;
    res.y_values[0] = y0;
    res.x_values[1] = x0 + h;
    res.y_values[1] = y1;

    BDF2Params params = { .h = h };

    for (int i = 2; i <= n; i++) {
        res.x_values[i] = res.x_values[i - 1] + h;

        params.y_n_prev = res.y_values[i - 2];
        params.y_n = res.y_values[i - 1];
        params.t_next = res.x_values[i];

        res.y_values[i] = secant(bdf2_equation, &params, params.y_n, params.y_n + 0.1, 1e-8, 1e-8);
    }

    return res;
}
