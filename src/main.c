#include <stdio.h>
#include <stdlib.h>
#include <methods.h>
#include "writeFile.h"

int main() {
    const double x0 = 0, y0 = 0;
    const double h = 0.25;
    const double n = 20;
    double y1 = 0.25 * f(0, 0);

    FILE *file = fopen("EDO_Descarga_de_Um_Capacitor_de_Circuitos_de_RC.txt", "w");

    if (file == NULL) {
        printf("Erro ao abrir arquivo para escrita!\n");
        return 1; 
    }

    fprintf(file, "Euler Explicito \n");
    Result res = euler(x0, y0, h, n);
    writeFile(file, res);  

    fprintf(file, "Euler Implicito \n");
    res = eulerImplicito(x0, y0, h, n);
    writeFile(file, res);  

    fprintf(file, "BDF2 \n");
    res = BDF2(x0, y1, y0, h, n);
    writeFile(file, res);  

    fprintf(file, "Adams Bashford 2 \n");
    res = adamsBashford2(x0, y0, h, n);
    writeFile(file, res);  

    fprintf(file, "Adams Moulton 2 \n");
    res = adamsMoulton2(y0, y1,x0, h, n);
    writeFile(file, res);  

    fprintf(file, "Trapezio Implicito \n");
    res = trapezioImplicito(x0, y0, h, n);
    writeFile(file, res);  

    fclose(file);
    free(res.x_values);
    free(res.y_values);

    return 0;
}