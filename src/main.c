#include <stdio.h>
#include <stdlib.h>
#include "methods.h"
#include "utils.h"
#include "error.h"
#include "ode_function.h"
#include "explicit_euler.h"
#include "implicit_euler.h"
#include "bdf2.h"
#include "adams_bashforth2.h"
#include "adams_moulton2.h"
#include "implicit_trapezoidal.h"

int main() {
    const double x0 = 0.0, y0 = 0.0;
    const double h = 0.25;
    const int n = 20;
    double y1 = y0 + h * ode_function(x0, y0);

    Result reference_solution = startEDO();

    FILE* file = fopen("RC_Circuit_Capacitor_Discharge.txt", "w");
    if (!file) {
        fprintf(stderr, "Error opening file for writing!\n");
        free(reference_solution.x_values);
        free(reference_solution.y_values);
        return 1;
    }

    Result res;

    fprintf(file, "Explicit Euler\n");
    res = explicit_euler(y0, x0, h, n);
    writeFile(file, res, error(res, reference_solution));
    free(res.x_values);
    free(res.y_values);

    fprintf(file, "Implicit Euler\n");
    res = implicit_euler(y0, x0, h, n);
    writeFile(file, res, error(res, reference_solution));
    free(res.x_values);
    free(res.y_values);

    fprintf(file, "BDF2\n");
    res = bdf2(y0, y1, x0, h, n);
    writeFile(file, res, error(res, reference_solution));
    free(res.x_values);
    free(res.y_values);

    fprintf(file, "Adams-Bashforth 2\n");
    res = adams_bashforth2(y0, x0, h, n);
    writeFile(file, res, error(res, reference_solution));
    free(res.x_values);
    free(res.y_values);

    fprintf(file, "Adams-Moulton 2\n");
    res = adams_moulton2(y0, y1, x0, h, n);
    writeFile(file, res, error(res, reference_solution));
    free(res.x_values);
    free(res.y_values);

    fprintf(file, "Implicit Trapezoidal\n");
    res = implicit_trapezoidal(y0, x0, h, n);
    writeFile(file, res, error(res, reference_solution));
    free(res.x_values);
    free(res.y_values);

    fprintf(file, "Reference Solution\n");
    writeFile(file, reference_solution, 0.0);

    fclose(file);
    free(reference_solution.x_values);
    free(reference_solution.y_values);

    return 0;
}