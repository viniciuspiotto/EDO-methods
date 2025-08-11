#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "error.h"
#include "methods.h"

double error (Result res, Result correct) {

    double soma = 0.0;
    double max_val = 0.0;

    for (int i = 0; i < correct.size; i++) {
        double val = fabs(correct.y_values[i]);
        if (val > max_val) {
            max_val = val;
        }
    }

    for (int i = 0; i < res.size; i++) {
        double diff = res.y_values[i] - correct.y_values[i];
        soma += pow(diff, 2);
    }

    double rms = sqrt(soma / res.size);

    if (max_val != 0.0) {
        return rms / max_val; 
    } else {
        return 0.0; 
    }
    
}
