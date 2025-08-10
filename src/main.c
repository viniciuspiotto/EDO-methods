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
    
    free(res.x_values);
    free(res.y_values);

    return 0;
}