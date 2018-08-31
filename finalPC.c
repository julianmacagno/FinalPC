#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define PI 3.141592653589793

void meshgrid(float* x, float* y, int tam_x, int tam_y, float** X, float** Y) {
    for(int i=0;i<tam_x; i++)
      X[i] = calloc(tam_x, sizeof(float));

    for(int i=0;i<tam_y; i++) 
      Y[i] = calloc(tam_y, sizeof(float));
    
    for(int i=0; i<tam_x; i++) {
        for(int j=0; j<tam_y; j++) {
            X[i][j] = x[j]; 
            //printf("%f  ", X[i][j]);
        }
        //printf("\n");
    }
    
    printf("\n");
    for(int i=0; i<tam_x; i++) {
        for(int j=0; j<tam_y; j++) {
            Y[i][j] = y[i];
            //printf("%f  ", Y[i][j]);
        }
        //printf("\n");
    }
}

void print(float *v, int tam);
void printM(float** v, int tam_x, int tam_y);
void cumtrapz();
float briereI(float tmps, float* pdes);
void varianza();
void vonFoerster(int nt, float* hmrs, float* tmps, float* pdes, float* t, float* tau);
void tempSim(float* t, int T0, int T1, int T365, float *tmeds, float *tmps, int tamVecs);
float* linspace(float start, float stop, int num);

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
    char fhrates[] = "BriereI";
    float* pdes = calloc(3, sizeof(float)); //np.zeros(3) 
    pdes[0] =  0.000131734;
    pdes[1] = 10.25740308;
    pdes[2] = 36.65400490;
    // Descripcion de la funcion de varianza en funcion de las temperaturas
    char fhnu[] = "Varianza";
    float* pnu = calloc(3, sizeof(float));//np.zeros(3) 
    pnu[0]  = 0.0;
    pnu[1]  = 0.0;
    pnu[2]  = 0.000223;
    // // Datos de la poblacion de ingreso Pulso 100 individuos 
    int* pinput  =  calloc(nt, sizeof(int));//np.zeros(nt)
    for(int i=0; i<4; i++) {
        pinput[i] = 25;
    }
    char idCorrida[] = "vFPy";
    int codfig = 8;
    vonFoerster(nt, hmrs, tmps, pdes, t, tau);
    // T, Tau, rates, RT, Pttau, ints, wts, pout = vonFoerster(dt, t, tau, nt, tmps, hmrs, pnu, fhnu, pdes, fhrates, pinput) 
    // //plt.plot(t, pinput, "b", t, pout, "r")
}

void cumtrapz() {

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

void varianza() {

}


void vonFoerster(int nt, float* hmrs, float* tmps, float* pdes, float* t, float* tau) {
    float tol = 0.0001;                      // von Foerster model tolerance
    print(t, nt);
    print(tau, nt);
    float** T = calloc(nt, sizeof(float *));
    float** Tau = calloc(nt, sizeof(float *));
    meshgrid(t, tau, nt, nt, T, Tau);          // T,Tau matrices //probablemente innecesaria, solo para graficar
    //printM(T, nt, nt); //vuelven con todos ceros, dentro de la funcion meshgrid esta todo bien
    //printM(Tau, nt, nt);
    //# Pttau   = np.zeros((nt,nt))           // initialize matrix which will hold the convolution kernel
    // wts     = np.zeros(nt)                // initialize vector which will hold weight for normalization
    float* rates = calloc(nt, sizeof(float));  // initialize vector of rates for each instant t
    float* hmrsnul = calloc(1, sizeof(float));
    for(int i=0; i<nt; i++) {
        if(abs(hmrs[i]) > 0)
        *hmrsnul += 1;
    }
    
    // evalua con humedad y sin humedad
    if (*hmrsnul != 0) {
       // calcula las tasas de desarrollo para temperaturas y humedades relativas dadas 
       for (int i=0; i<nt; i++) {
            //rates[i] = briereI(tmps[i], hmrs[i], pdes);   //lamar a BRIEREI //no se usa nunca
       }
    }
    else {
        // calcula las tasas de desarrollo  para temperaturas dadas
        for (int i=0; i<nt; i++) {
            rates[i] = briereI(tmps[i], pdes); //lamar a BRIEREI
        }
    }
             
//     RT    = np.dot(dt * np.ones((nt,1)), [cumtrapz(rates)])  // create a matrix which is the cumulative development
//     RTau  = np.transpose(RT)                                    // create transpose of cumulative matrix for use in kernel
//     vexp  = np.vectorize(exp)   
//     vsqrt = np.vectorize(sqrt)    
    
//     if (pnu[0] == 0) and (pnu[1] == 0):
//         nu = pnu[2]
//         Pttau  = (T>Tau)*vexp(-(1-(RT-RTau))**2/(4*nu*(abs(T-Tau)+tol)))\
//                            /vsqrt(4*pi*nu*(abs(T-Tau)**3+tol))         // extended von foerster kernel
//     else:
//         nus   = feval(fhnu, tmps, pnu)      // calcula las varianzas para temperatures dadas //lamar a varianza
//         NU    = np.dot(np.ones(nt,1), nus)// crea una matriz de varianzas en funcion de las temperaturas
//         Pttau = (T>Tau)*exp(-(1-(RT-RTau))**2/(4*NU*(abs(T-Tau)+tol))) \
//                            /vsqrt(4*pi*NU*(abs(T-Tau)**3+tol))         // extended von foerster kernel
    
//     ints = dt*np.trapz(Pttau,axis=1)    // integrate in columns to normalize
//     wts  = (ints>tol)*ints+(ints<=tol)  // calculate a  weighting factor 
//                                          // make it one if the integral is too  small (< tol)
//     pout = dt*np.dot(pinput/np.transpose(wts), Pttau)      // output distribution, normalized by integral of P
}

//TEMPSIM simula tabla de Temperaturas para hemisferio sur
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
            printf("%f  ", &v[1][1]);
        }
        printf("\n");
    }
    printf("\n\n");
}