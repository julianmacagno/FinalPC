#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define PI 3.141592653589793

void print(float *v, int tam);
void printM(float** v, int tam_x, int tam_y);
void cumtrapz(float *v1, float *v2, int tam);
float briereI(float tmps, float* pdes);
void vonFoerster(int nt, float* hmrs, float* tmps, float* pdes, float* t, float* tau, float dt);
void tempSim(float* t, int T0, int T1, int T365, float *tmeds, float *tmps, int tamVecs);
float* linspace(float start, float stop, int num);
float** initZeros2DFloatMatrix(int tam_x, int tam_y) {
    float **ptr = calloc(tam_x, sizeof(float*));
    for(int i=0;i<tam_y; i++)
      ptr[i] = calloc(tam_y, sizeof(float));
    
    return ptr;
}

void meshgrid(float* x, float* y, int tam_x, int tam_y, float** X, float** Y) {
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

void main() {
    int tmin = 30;         // dia inicial de la corrida 30 de enero
    int tmax = 30 + 100;   // dia final de la corrida 30 de enero mas 100 dias
    int nd = tmax - tmin;  // numero de dias
    int td = 4;            // intervalos de tiempo diarios
    int nt = nd * td;      // numero total de pasos de tiempo
    float dt = (float)nd / nt; // paso discreto de tiempo
    float *t = linspace(tmin,tmax-dt,nt); // vector de tiempos, t devuelve un vector de numeros entre un rango separados por el rango dividido la cantidad de intervalos
    float *tau  = t;            // vector of tiempos, tau
    int T0   = 15;         // T0   : temperatura media
    int T1   = 15;         // T1   : amplitud termica anual
    int T365 = 15;         // T365 : amplitud termica diaria 
    float *tmeds = calloc(nt, sizeof(float)); //se crean aca los vectores del tamaño correspondiente
    float* tmps = calloc(nt, sizeof(float));
    tempSim(t, T0, T1, T365, tmeds, tmps, nt);
    float* hmrs  = calloc(nt, sizeof(float)); //reemplazo de //hmrs = np.zeros(nt)
    // Datos para huevos de Diatrea 
    float* pdes = calloc(3, sizeof(float)); //np.zeros(3) 
    pdes[0] =  0.000131734;
    pdes[1] = 10.25740308;
    pdes[2] = 36.65400490;
    // Descripcion de la funcion de varianza en funcion de las temperaturas
    float* pnu = calloc(3, sizeof(float));//np.zeros(3) 
    pnu[0]  = 0.0;
    pnu[1]  = 0.0;
    pnu[2]  = 0.000223;
    // // Datos de la poblacion de ingreso Pulso 100 individuos 
    int* pinput  =  calloc(nt, sizeof(int));//np.zeros(nt)
    for(int i=0; i<4; i++) {
        pinput[i] = 25;
    }
    int codfig = 8;
    vonFoerster(nt, hmrs, tmps, pdes, t, tau, dt);
    // T, Tau, rates, RT, Pttau, ints, wts, pout = vonFoerster(dt, t, tau, nt, tmps, hmrs, pnu, fhnu, pdes, fhrates, pinput) 
    // //plt.plot(t, pinput, "b", t, pout, "r")
}

float briereI(float tmpsi, float* p) {
    float r  = 0;
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

void vonFoerster(int nt, float* hmrs, float* tmps, float* pdes, float* t, float* tau, float dt) {
    float tol = 0.0001;                         //von Foerster model tolerance
    
    float** T = initZeros2DFloatMatrix(nt, nt);             //Create T pointer to 2x2 matrix
    float** Tau = initZeros2DFloatMatrix(nt, nt);           //Create T pointer to 2x2 matrix
    meshgrid(t, tau, nt, nt, T, Tau);           // T,Tau matrices
    
    float** Pttau = initZeros2DFloatMatrix(nt, nt);         // initialize matrix which will hold the convolution kernel
    float* wts = calloc(nt, sizeof(float));     // initialize vector which will hold weight for normalization
    float* rates = calloc(nt, sizeof(float));   // initialize vector of rates for each instant t
    
    float* hmrsnul = calloc(1, sizeof(float));
    for(int i=0; i<nt; i++) {
        if(abs(hmrs[i]) > 0)
        *hmrsnul += 1;
    }

    // calcula las tasas de desarrollo  para temperaturas dadas //EX IF-ELSE
    for (int i=0; i<nt; i++) {
        rates[i] = briereI(tmps[i], pdes); //lamar a BRIEREI
    }

    float *vec = calloc(nt, sizeof(float));
    cumtrapz(rates, vec, nt); //returns vec
      
    float **RT = initZeros2DFloatMatrix(nt, nt);
    for(int i=0; i<nt; i++) {
        for(int j=0; j<nt; j++) {
            RT[i][j] = dt * vec[j];
        }
    }
//     RTau  = np.transpose(RT)                                    // create transpose of cumulative matrix for use in kernel
//     vexp  = np.vectorize(exp)   
//     vsqrt = np.vectorize(sqrt)    
    
//     nu = pnu[2] //EX IF-ELSE
//     Pttau  = (T>Tau)*vexp(-(1-(RT-RTau))**2/(4*nu*(abs(T-Tau)+tol)))\
//                  /vsqrt(4*pi*nu*(abs(T-Tau)**3+tol))         // extended von foerster kernel
    
//     ints = dt*np.trapz(Pttau,axis=1)    // integrate in columns to normalize
//     wts  = (ints>tol)*ints+(ints<=tol)  // calculate a  weighting factor 
//                                          // make it one if the integral is too  small (< tol)
//     pout = dt*np.dot(pinput/np.transpose(wts), Pttau)      // output distribution, normalized by integral of P
}

void cumtrapz(float *v1, float *v2, int tam) {
    v2[0]=0;
    for(int i=1; i<tam; i++) {
        v2[i]= v2[i-1] + (v1[i-1] + v1[i]) / 2;
    }
}

void tempSim(float* t, int T0, int T1, int T365, float* tmeds, float* tmps, int tamVecs) {    
    for(int i=0; i<tamVecs; i++) {
        tmeds[i] = T0 + T1/2 * cos(2 * (PI/365) * t[i]); //equivalente a hacer np.factorize y desp resolver con el vector t. Hay un problema con la precision de los calculos, en python tiene mas decimales
    }

    for(int i=0; i<tamVecs; i++) {
        tmps[i] = tmeds[i] - T365/2 * cos(2 * PI * (t[i] - (int)t[i])); //equivalente a hacer np.factorize y desp resolver con el vector t y tmeds. Hay un problema con la precision de los calculos, en python tiene mas decimales
    }
}

float* linspace(float start, float stop, int num) {
    float *v = calloc(num, sizeof(float));
    float val = (stop-start)/(num-1); //Agregue el -1 para que los resultados sean los mismos que en el programa de python
    for(int i=0; i<num; i++) {
        v[i] = start + i*val;
    }

    return v;
}

void print(float *v, int tam) {
    for(int i=0; i<tam; i++) 
        printf("%f  \n", v[i]);
}

void printM(float** v, int tam_x, int tam_y) {
    for(int i=0; i<tam_x; i++) {
        for(int j=0; j<tam_y; j++) {
            printf("%f  ", v[i][j]);
        }
        printf("\n");
    }
    printf("\n\n");
}
