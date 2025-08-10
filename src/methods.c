#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <methods.h>

double g_y_n;
double g_t_next;
double g_h;

double f (double x, double y){
    return 10 * (12 * cos(2 * M_PI * x) - y);
}

double g (double y) {
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

double secante(double x0, double x1, double epsilon, double delta) {

    double x2 = x1;  
    int zeroDivisionError = 0;

    while (fabs(x2 - x1) > delta || fabs(g(x2)) > epsilon) {
        if (fabs(g(x1) - g(x0)) == 0) {
            zeroDivisionError = 1;
            break;
        } else {
            x2 = x1 - g(x1) * (x1 - x0) / (g(x1) - g(x0));
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

        res.y_values[i] = secante(g_y_n, g_y_n + 0.1, 1e-8, 1e-8);
    }

    return res;

}

