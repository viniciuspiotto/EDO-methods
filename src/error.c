#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "error.h"
#include "methods.h"

double iterationError(Result res, Result correct, int start, int end){

    double soma = 0.0;
    double max_val = 0.0;

    for (int i = start; i < end; i++) {
        double val = fabs(correct.y_values[i]);
        if (val > max_val) {
            max_val = val;
        }
    }

    for (int i = start; i < end; i++) {
        double diff = res.y_values[i] - correct.y_values[i];
        soma += pow(diff, 2);
    }

    double rms = sqrt(soma / (end - start));

    if (max_val != 0.0) {
        return rms / max_val; 
    } else {
        return 0.0; 
    }

}

double * error(Result res, Result correct) {

    double * errors = malloc(5 * sizeof(double));

    for (int i = 0; i < 4; i++){
       errors[i] = iterationError(res, correct, i * 5, (i != 3) ? (i + 1) * 5 : ((i + 1) * 5) + 1);
    }

    errors[4] = iterationError(res, correct, 0, 21);
    return errors;
    
}
