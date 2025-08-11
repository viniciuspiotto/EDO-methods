#include <stdio.h>
#include <stdlib.h>
#include <methods.h>

int main() {

    printf("euler explicito:\n");
    Result res = euler(0, 0, 0.25, 20);
    for(int i = 0; i < res.size; i++){
        printf("x: %.4f y: %.4f\n", res.x_values[i], res.y_values[i]);
    }

    printf("euler implicito:\n");
    res = eulerImplicito(0, 0, 0.25, 20);
    for(int i = 0; i < res.size; i++){
        printf("x: %.4f y: %.4f\n", res.x_values[i], res.y_values[i]);
    }

    printf("BDF2:\n");
    res = BDF2(0, 0.25 * f(0, 0), 0, 0.25, 20);
    for(int i = 0; i < res.size; i++){
        printf("x: %.4f y: %.4f\n", res.x_values[i], res.y_values[i]);
    }

    printf("Adam Bashford 2: \n");
    res = adamsBashford2(0, 0, 0.25, 20);
    for(int i = 0; i < res.size; i++){
        printf("x: %.4f y: %.4f\n", res.x_values[i], res.y_values[i]);
    }
    
    free(res.x_values);
    free(res.y_values);

    return 0;
}