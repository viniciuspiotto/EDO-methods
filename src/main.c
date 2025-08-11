#include <stdio.h>
#include <stdlib.h>
#include <methods.h>
#include "writeFile.h"

int main() {

    FILE *file = fopen("EDO_Descarga_de_Um_Capacitor_de_Circuitos_de_RC.txt", "w");

    if (file == NULL) {
        printf("Erro ao abrir arquivo para escrita!\n");
        return 1; 
    }

    fprintf(file, "Euler Explicito \n");
    Result res = euler(0, 0, 0.25, 20);
    writeFile(file, res);  

    fprintf(file, "Euler Implicito \n");
    res = eulerImplicito(0, 0, 0.25, 20);
    writeFile(file, res);  

    fprintf(file, "BDF2 \n");
    res = BDF2(0, 0.25 * f(0, 0), 0, 0.25, 20);
    writeFile(file, res);  

    fprintf(file, "Adams Bashford 2 \n");
    res = adamsBashford2(0, 0, 0.25, 20);
    writeFile(file, res);  

    fprintf(file, "Trapezio Implicito \n");
    res = trapezioImplicito(0, 0, 0.25, 20);
    writeFile(file, res);  

    fclose(file);
    free(res.x_values);
    free(res.y_values);

    return 0;
}