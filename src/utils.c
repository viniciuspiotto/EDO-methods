#include <stdio.h>
#include <stdlib.h>
#include <methods.h>
#include "utils.h"

Result expected_result() {
    Result correct = {.size = 21};

    correct.x_values = malloc(correct.size * sizeof(double));
    correct.y_values = malloc(correct.size * sizeof(double));
    if (!correct.x_values || !correct.y_values) {
        fprintf(stderr, "Memory allocation failed in expected_result\n");
        exit(EXIT_FAILURE);
    }

    double x_tmp[21] = {
        0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25,
        2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0
    };

    double y_tmp[21] = {
        0, 4.69951014480389, -8.66145140695345, -5.41048537162463, 8.60309100643824,
        5.40569485822865, -8.60348423572376, -5.40572713645401, 8.60348158616568,
        5.40572691896504, -8.60348160401826, -5.40572692043047, 8.60348160389797,
        5.40572692042060, -8.60348160389878, -5.40572692042066, 8.60348160389877,
        5.40572692042066, -8.60348160389877, -5.40572692042066, 8.60348160389877
    };

    for (int i = 0; i < correct.size; i++) {
        correct.x_values[i] = x_tmp[i];
        correct.y_values[i] = y_tmp[i];
    }

    return correct;
}

void print_result(Result res) {
    for (int i = 0; i < res.size; i++) {
        printf("    %.10f %.10f\n", res.x_values[i], res.y_values[i]);
    }
}

void write_result(FILE* file, Result res) {
    for (int i = 0; i < res.size; i++) {
        fprintf(file, "    %.10f %.10f\n", res.x_values[i], res.y_values[i]);
    }
}