#include "adams_moulton2.h"
#include "ode_function.h"
#include "secant.h"
#include <stdlib.h>

static double adams_moulton2_equation(double y, void* params) {
    AdamsMoulton2Params* p = (AdamsMoulton2Params*) params;
    return 12 * (y - p->y_n) - p->h * (5 * ode_function(p->t_next, y) + 8 * p->f_k - p->f_k_1);
}

Result adams_moulton2(double y0, double y1, double x0, double h, int n) {
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

    double fk = ode_function(res.x_values[0], res.y_values[0]);
    double fk_1 = ode_function(res.x_values[1], res.y_values[1]);

    AdamsMoulton2Params params = {
        .h = h,
        .y_n = res.y_values[0],
        .t_next = res.x_values[1],
        .f_k = fk_1,
        .f_k_1 = fk
    };

    for (int i = 2; i <= n; i++) {
        res.x_values[i] = res.x_values[i - 1] + h;

        params.h = h;
        params.t_next = res.x_values[i];
        params.y_n = res.y_values[i - 1];
        params.f_k = fk_1;
        params.f_k_1 = fk;

        double initial_guess = res.y_values[i - 1];
        res.y_values[i] = secant(adams_moulton2_equation, &params, initial_guess, initial_guess + 0.1, 1e-8, 1e-8);

        fk = fk_1;
        fk_1 = ode_function(res.x_values[i], res.y_values[i]);
    }

    return res;
}