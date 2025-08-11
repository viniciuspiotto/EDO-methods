#include <stdio.h>
#include <stdlib.h>
#include <methods.h>
#include "writeFile.h"

void writeFile (FILE * file, Result res){

    if (file == NULL) return;  
    
    for (int i = 0; i < res.size; i++) {
        fprintf(file, "%.10f %.10f\n", res.x_values[i], res.y_values[i]);
    }

    fprintf(file, "\n");

}