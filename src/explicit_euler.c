#include "explicit_euler.h"
#include "ode_function.h"
#include <stdlib.h>

Result explicit_euler(double y0, double x0, double h, int n) {
    int size = n + 1;
    Result res = {
        .size = size,
        .x_values = malloc(size * sizeof(double)),
        .y_values = malloc(size * sizeof(double))
    };

    res.x_values[0] = x0;
    res.y_values[0] = y0;

    for (int i = 1; i <= n; i++) {
        res.x_values[i] = res.x_values[i - 1] + h;
        res.y_values[i] = res.y_values[i - 1] + h * ode_function(res.x_values[i - 1], res.y_values[i - 1]);
    }

    return res;
}