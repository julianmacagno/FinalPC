#include <stdio.h>
#include <stdlib.h>

void meshgrid(float* x, float* y, int tam_x, int tam_y, float** X, float** Y) {
    X = calloc(tam_x, sizeof(float *));
    for(int i=0;i<tam_x; i++)
      X[i] = calloc(tam_x, sizeof(float));

    Y = calloc(tam_y, sizeof(float *));
    for(int i=0;i<tam_y; i++) 
      Y[i] = calloc(tam_y, sizeof(float));
    
    for(int i=0; i<tam_x; i++) {
        for(int j=0; j<tam_y; j++) {
            X[i][j] = x[j]; 
            printf("%f  ", X[i][j]);
        }
        printf("\n");
    }
    
    printf("\n");
    for(int i=0; i<tam_x; i++) {
        for(int j=0; j<tam_y; j++) {
            Y[i][j] = y[i];
            printf("%f  ", Y[i][j]);
        }
        printf("\n");
    }
}

float* linspace(float start, float stop, int num) {
    float *v = calloc(num, sizeof(float));
    float val = (stop-start)/(num-1); //Agregue el -1 para que los resultados sean los mismos que en el programa de python
    for(int i=0; i<num; i++) {
        v[i] = start + i*val;
        printf("%f  ", v[i]);
    }
    printf("\n");
    return v;
}

int main() {
    float *v = calloc(10, sizeof(float));
    v = linspace(1, 10, 10);
    printf("\n");
    float *b = calloc(10, sizeof(float));
    b = linspace(1, 10, 10);
    printf("\n");
    float** X, **Y;
    meshgrid(v, b, 10, 10, X, Y);
    printf("\n");
}