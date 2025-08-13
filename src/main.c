#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
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
#include "runge_kutta2.h"

int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <method>\n", argv[0]);
        fprintf(stderr, "Methods:\n");
        fprintf(stderr, "1 - Explicit Euler\n");
        fprintf(stderr, "2 - Implicit Euler\n");
        fprintf(stderr, "3 - BDF2\n");
        fprintf(stderr, "4 - Adams-Bashforth 2\n");
        fprintf(stderr, "5 - Adams-Moulton 2\n");
        fprintf(stderr, "6 - Implicit Trapezoidal\n");
        fprintf(stderr, "7 - Runge Kutta 2\n");
        return 1;
    }

    int choice = atoi(argv[1]);
    if (choice < 1 || choice > 7) {
        fprintf(stderr, "Invalid option: %d\n", choice);
        return 1;
    }

    const double x0 = 0.0, y0 = 0.0;
    const double h = 0.25;
    const int n = 20;
    double y1 = y0 + h * ode_function(x0, y0);
    Result reference_solution = expected_result();

    Result res;
    const char *method_name = NULL;

    struct timespec t_start, t_end;
    clock_gettime(CLOCK_MONOTONIC_RAW, &t_start);

    switch (choice) {
        case 1: 
            method_name = "Explicit Euler";
            res = explicit_euler(y0, x0, h, n);
            break;
        case 2: 
            method_name = "Implicit Euler";
            res = implicit_euler(y0, x0, h, n);
            break;
        case 3: 
            method_name = "BDF2";
            res = bdf2(y0, y1, x0, h, n);
            break;
        case 4: 
            method_name = "Adams-Bashforth 2";
            res = adams_bashforth2(y0, x0, h, n);
            break;
        case 5: 
            method_name = "Adams-Moulton 2";
            res = adams_moulton2(y0, y1, x0, h, n);
            break;
        case 6: 
            method_name = "Implicit Trapezoidal";
            res = implicit_trapezoidal(y0, x0, h, n);
            break;
        case 7: 
            method_name = "Runge Kutta 2";
            res = runge_kutta2(y0, x0, h, n);
            break;
    }

    clock_gettime(CLOCK_MONOTONIC_RAW, &t_end);

    double elapsed_time = (t_end.tv_sec - t_start.tv_sec) + (t_end.tv_nsec - t_start.tv_nsec) / 1e9;
    double * err = error(res, reference_solution);

    printf("%s\n", method_name);
    printf("Elapsed Time: %.15f s\n", elapsed_time);
    printf("Error: %.15f\n", err[4]);
    for(int i = 0; i < 4; i++){
        printf("Error Q%d: %.15f\n", i+1, err[i]);
    }
    printf("Results:\n");
    print_result(res);

    FILE *file = fopen("output.txt", "a");
    if (!file) {
        fprintf(stderr, "Error opening output.txt for writing!\n");
        free(reference_solution.x_values);
        free(reference_solution.y_values);
        free(res.x_values);
        free(res.y_values);
        return 1;
    }
    fprintf(file, "%s\n", method_name);
    fprintf(file, "Elapsed Time: %.15f s\n", elapsed_time);
    fprintf(file, "Error: %.15f\n", err[4]);
    for(int i = 0; i < 4; i++){
        fprintf(file, "Error Q%d: %.15f\n", i+1, err[i]);
    }
    fprintf(file, "Results:\n");
    write_result(file, res);
    fprintf(file, "\n");

    fclose(file);

    free(reference_solution.x_values);
    free(reference_solution.y_values);
    free(res.x_values);
    free(res.y_values);

    return 0;
}