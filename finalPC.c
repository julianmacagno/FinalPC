#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#define PI 3.141592653589793
#define e 2.71828182846

void cumtrapz(double *v1, double *v2, int tam);
double briereI(double tmps, double* pdes);
void vonFoerster(int nt, double* hmrs, double* tmps, double* pdes, double* t, double* tau, double dt, double* pnu, double* pinput);
void tempSim(double* t, int T0, int T1, int T365, double *tmeds, double *tmps, int tamVecs);
double* linspace(double start, double stop, int num);
void findTranspose(double **vec1, double** vec2, int x, int y);
double** initZeros2DFloatMatrix(int tam_x, int tam_y);
double valAbs(double a) { return (a<0) ? a*-1 : a; }
double* trapz(double** mat, int tam);
void meshgrid(double* x, double* y, int tam_x, int tam_y, double** X, double** Y);
void printV(double *v, int tam);
void printM(double** v, int tam_x, int tam_y);

void main() {
    clock_t begin = clock();
    int tmin = 30;         // dia inicial de la corrida 30 de enero
    int tmax = 30 + 100;   // dia final de la corrida 30 de enero mas 100 dias
    int nd = tmax - tmin;  // numero de dias
    int td = 4;            // intervalos de tiempo diarios
    int nt = nd * td;      // numero total de pasos de tiempo
    double dt = (double)nd / nt; // paso discreto de tiempo
    double *t = linspace(tmin,tmax-dt,nt); // vector de tiempos, t devuelve un vector de numeros entre un rango separados por el rango dividido la cantidad de intervalos
    double *tau  = t;            // vector of tiempos, tau
    int T0   = 15;         // T0   : temperatura media
    int T1   = 15;         // T1   : amplitud termica anual
    int T365 = 15;         // T365 : amplitud termica diaria 
    double *tmeds = calloc(nt, sizeof(double)); //se crean aca los vectores del tamaño correspondiente
    double* tmps = calloc(nt, sizeof(double));
    tempSim(t, T0, T1, T365, tmeds, tmps, nt);
    double* hmrs  = calloc(nt, sizeof(double)); //reemplazo de //hmrs = np.zeros(nt)
    // Datos para huevos de Diatrea 
    double* pdes = calloc(3, sizeof(double)); //np.zeros(3) 
    pdes[0] =  0.000131734;
    pdes[1] = 10.25740308;
    pdes[2] = 36.65400490;
    // Descripcion de la funcion de varianza en funcion de las temperaturas
    double* pnu = calloc(3, sizeof(double));//np.zeros(3) 
    pnu[0]  = 0.0;
    pnu[1]  = 0.0;
    pnu[2]  = 0.000223;
    // // Datos de la poblacion de ingreso Pulso 100 individuos 
    double* pinput  =  calloc(nt, sizeof(int));//np.zeros(nt)
    for(int i=0; i<4; i++) {
        pinput[i] = 25;
    }
    int codfig = 8;
    vonFoerster(nt, hmrs, tmps, pdes, t, tau, dt, pnu, pinput);
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Time elapsed: %g miliseconds\n", time_spent*1000);
}

void vonFoerster(int nt, double* hmrs, double* tmps, double* pdes, double* t, double* tau, double dt, double* pnu, double* pinput) {
    double tol = 0.0001;                         //von Foerster model tolerance
    
    double** T = initZeros2DFloatMatrix(nt, nt);             //Create T pointer to 2x2 matrix
    double** Tau = initZeros2DFloatMatrix(nt, nt);           //Create T pointer to 2x2 matrix
    meshgrid(t, tau, nt, nt, T, Tau);           // T,Tau matrices
    
    double* rates = calloc(nt, sizeof(double));   // initialize vector of rates for each instant t
    
    double* hmrsnul = calloc(1, sizeof(double));
    for(int i=0; i<nt; i++) {
        if(abs(hmrs[i]) > 0)
        *hmrsnul += 1;
    }

    for (int i=0; i<nt; i++) { // calcula las tasas de desarrollo  para temperaturas dadas
        rates[i] = briereI(tmps[i], pdes); //lamar a BRIEREI
    }

    double *vec = calloc(nt, sizeof(double));
    cumtrapz(rates, vec, nt); //returns vec

    double **RT = initZeros2DFloatMatrix(nt, nt);
    for(int i=0; i<nt; i++) {
        for(int j=0; j<nt; j++) {
            RT[i][j] = dt * vec[j];
        }
    }

    double** RTau = initZeros2DFloatMatrix(nt, nt);
    findTranspose(RT, RTau, nt, nt); // create transpose of cumulative matrix for use in kernel
    double nu = pnu[2];
    double** Pttau = initZeros2DFloatMatrix(nt, nt);         // initialize matrix which will hold the convolution kernel
    for(int i=0; i<nt; i++) {
        for(int j=0; j<nt; j++) {
            double a = (T[i][j] > Tau[i][j]);
            double b = - pow(1 - (RT[i][j] - RTau[i][j]), 2);
            double c = (4 * nu * (valAbs(T[i][j] - Tau[i][j]) + tol));
            double numerador =  pow(e, b / c);
            double denominador =  sqrt(4 * PI * nu * (pow(valAbs(T[i][j] - Tau[i][j]), 3) + tol));
            Pttau[i][j] = a * numerador / denominador;
        }
    }
    double* ints = trapz(Pttau, nt);
    for(int i=0; i<nt; i++) {
        ints[i] *= dt;
    }

    double* wts = calloc(nt, sizeof(double));     // initialize vector which will hold weight for normalization
    for(int i=0; i<nt; i++) {
        wts[i] = (ints[i]>tol)*ints[i] + (ints[i]<=tol); // calculate a  weighting factor 
    }
    
    double* pout = calloc(nt, sizeof(double));
    for (int i = 0; i < nt; i++) {
        double result = 0.0;
        for (int j = 0; j < nt; j++) {
            result += Pttau[j][i] * (pinput[j]/wts[j]);
        }
        pout[i] = dt * result;
    }
}

void findTranspose(double **vec, double** transpose, int x, int y) { //return transpose of a matrix on **transpose
    for(int i=0; i<x; ++i) {
        for(int j=0; j<y; ++j) {
            transpose[j][i] = vec[i][j];
        }
    }
}

double briereI(double tmpsi, double* p) {
    double r  = 0;
	if (tmpsi <= p[1]) {
	   r = 0;
    }
	else {
        if ((tmpsi > p[1]) && (tmpsi <= p[2])) {
            r = p[0] * tmpsi * (tmpsi - p[1]) * sqrt(p[2] - tmpsi);
        }	
        else {
           r = 0;
        }
    }
	return r;
}

void cumtrapz(double *v1, double *v2, int tam) {
    v2[0]=0;
    for(int i=1; i<tam; i++) {
        v2[i]= v2[i-1] + (v1[i-1] + v1[i]) / 2;
    }
}

void tempSim(double* t, int T0, int T1, int T365, double* tmeds, double* tmps, int tamVecs) {    
    for(int i=0; i<tamVecs; i++) {
        tmeds[i] = T0 + T1/2 * cos(2 * (PI/365) * t[i]); //equivalente a hacer np.factorize y desp resolver con el vector t. Hay un problema con la precision de los calculos, en python tiene mas decimales
    }

    for(int i=0; i<tamVecs; i++) {
        tmps[i] = tmeds[i] - T365/2 * cos(2 * PI * (t[i] - (int)t[i])); //equivalente a hacer np.factorize y desp resolver con el vector t y tmeds. Hay un problema con la precision de los calculos, en python tiene mas decimales
    }
}

double* linspace(double start, double stop, int num) {
    double *v = calloc(num, sizeof(double));
    double val = (stop-start)/(num-1); //Agregue el -1 para que los resultados sean los mismos que en el programa de python
    for(int i=0; i<num; i++) {
        v[i] = start + i*val;
    }

    return v;
}

void printV(double *v, int tam) {
    for(int i=0; i<tam; i++) {
        printf("%g  \n", v[i]);
    }
}

void printM(double** v, int tam_x, int tam_y) {
    for(int i=0; i<tam_x; i++) {
        for(int j=0; j<tam_y; j++) {
            printf("%g  ", v[i][j]);
        }
        printf("\n");
    }
    printf("\n\n");
}

double** initZeros2DFloatMatrix(int tam_x, int tam_y) {
    double **ptr = calloc(tam_x, sizeof(double*));
    for(int i=0;i<tam_y; i++) {
      ptr[i] = calloc(tam_y, sizeof(double));
    }
    
    return ptr;
}

double* trapz(double** mat, int tam) {
    double *res=calloc(tam, sizeof(double));
    double sum = 0;
    for(int i=0; i<tam; i++) {
        for(int j=0; j<tam-1; j++) {
           sum += (mat[i][j+1] + mat[i][j]) / 2.0;
        }
        res[i] = sum;
        sum=0;
    }
    return res;    
}

void meshgrid(double* x, double* y, int tam_x, int tam_y, double** X, double** Y) {
    for(int i=0; i<tam_x; i++) {
        for(int j=0; j<tam_y; j++) {
            X[i][j] = x[j]; 
        }
    }

    for(int i=0; i<tam_x; i++) {
        for(int j=0; j<tam_y; j++) {
            Y[i][j] = y[i];
        }
    }
}