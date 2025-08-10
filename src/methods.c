#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <methods.h>

double f (double x, double y){
    return 10 * (12 * cos(2 * PI * x) - y);
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

// def adams_bashforth_2(f, y0, t0, h, n):
//     y = [y0]
//     t = [t0]
//     y.append(y0 + h * f(t0, y0))
//     t.append(t0 + h)

//     for k in range(1, n):
//         fk = f(t[k], y[k])
//         fk_1 = f(t[k-1], y[k-1])
//         y_next = y[k] + h * (3/2 * fk - 1/2 * fk_1)
//         y.append(y_next)
//         t.append(t[k] + h)
//     return t, y