#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define PI 3.141592653589793

void print(float *v, int tam) {
    for(int i=0; i<tam; i++) 
        printf("%f  ", v[i]);
}
void feval();
void cumtrapz();
void briereI();
void tempSim(float* t, int T0, int T1, int T365, float *tmeds, float *tmps, int tamVecs);
void varianza();
void vonFoerster();
float* linspace(float start, float stop, int num);

void main() {
    //Generacion de datos de prueba de una corrida 
    int tmin = 30;         // dia inicial de la corrida 30 de enero
    int tmax = 30 + 100;   // dia final de la corrida 30 de enero mas 100 dias
    int nd = tmax - tmin;  // numero de dias
    int td = 4;            // intervalos de tiempo diarios
    int nt = nd * td;      // numero total de pasos de tiempo
    float dt = (float)nd / nt; // paso discreto de tiempo
    float *t = linspace(tmin,tmax-dt,nt); // vector de tiempos, t devuelve un vector de numeros entre un rango separados por el rango dividido la cantidad de intervalos
    // tau  = t            // vector of tiempos, tau
    int T0   = 15;         // T0   : temperatura media
    int T1   = 15;         // T1   : amplitud termica anual
    int T365 = 15;         // T365 : amplitud termica diaria 
    float *tmeds, *tmps; //necesarios para devolver el valor
    tempSim(t, T0, T1, T365, tmeds, tmps, nt);
    float* hmrs  = calloc(nt, sizeof(float)); //reemplazo de #hmrs = np.zeros(nt)
    // Datos para huevos de Diatrea 
    // fhrates = "BriereI"
    // pdes    = np.zeros(3) 
    // pdes[0] =  0.000131734
    // pdes[1] = 10.25740308  
    // pdes[2] = 36.65400490
    // // Descripcion de la funcion de varianza en funcion de las temperaturas
    // fhnu    = "Varianza"
    // pnu     = np.zeros(3) 
    // pnu[0]  = 0.0
    // pnu[1]  = 0.0  
    // pnu[2]  = 0.000223
    // // Datos de la poblacion de ingreso Pulso 100 individuos 
    // pinput  =  np.zeros(nt)
    // for i in range(4):
    // pinput[i] = 25
    // idCorrida = "vFPy"
    // codfig = 8
    // T, Tau, rates, RT, Pttau, ints, wts, pout = vonFoerster(dt, t, tau, nt, tmps, hmrs, pnu, fhnu, pdes, fhrates, pinput) 
    // //plt.plot(t, pinput, "b", t, pout, "r")
}

void feval() {}
void cumtrapz() {}
void briereI() {}
void varianza() {}
void vonFoerster() {}

//TEMPSIM simula tabla de Temperaturas para hemisferio sur
void tempSim(float* t, int T0, int T1, int T365, float* tmeds, float* tmps, int tamVecs) {    
    tmeds = calloc(tamVecs, sizeof(float)); //se crean aca los vectores del tamaño correspondiente
    tmps = calloc(tamVecs, sizeof(float));
    for(int i=0; i<tamVecs; i++) {
        tmeds[i] = T0 + T1/2 * cos(2 * (PI/365) * t[i]); //equivalente a hacer np.factorize y desp resolver con el vector t. Hay un problema con la precision de los calculos, en python tiene mas decimales
    }

    for(int i=0; i<tamVecs; i++) {
        tmps[i] = tmeds[i] - T365/2 * cos(2 * PI * (t[i] - (int)t[i])); //equivalente a hacer np.factorize y desp resolver con el vector t y tmeds. Hay un problema con la precision de los calculos, en python tiene mas decimales
    }

    // for(int i=0; i<tamVecs; i+=6)
    //     printf("%f %f %f %f %f %f\n", tmps[i], tmps[i+1], tmps[i+2], tmps[i+3], tmps[i+4], tmps[i+5]);
}

float* linspace(float start, float stop, int num) {
    float *v = calloc(num, sizeof(float));
    float val = (stop-start)/(num-1); //Agregue el -1 para que los resultados sean los mismos que en el programa de python
    for(int i=0; i<num; i++) {
        v[i] = start + i*val;
    }

    return v;
}