import numpy as np
from math import *
import timeit

def feval(funcName  , *args):      # first argument must be a string
    return eval(funcName)(*args)

def cumtrapz(hs):
    ct = np.zeros(len(hs))
    for i in range(1,len(hs)):
        ct[i] = ct[i-1] + (hs[i-1] + hs[i]) / 2.0
    return ct
# BRIEREI - BriereI - 3 parametros      
def BriereI(tmpsi, p):
	r  = 0.0
	if tmpsi <= p[1]:
	   r = 0.0
	elif (tmpsi > p[1]) and (tmpsi <= p[2]):
	   r = p[0] * tmpsi * (tmpsi - p[1]) * sqrt(p[2] - tmpsi)
	else:
	   r = 0.0
	return r
	
def  TempSim(t, T0, T1, T365):
#TEMPSIM simula tabla de Temperaturas para hemisferio sur
# t    : vector de tiempos en dias.fraccion desde t[0] instante inicial
# T0   : temperatura media
# T1   : amplitud termica anual
# T365 : amplitud termica diaria
    vcos  = np.vectorize(cos)
    vint  = np.vectorize(int)
    T1    = T1/2
    T365  = T365/2
    tmeds = T0+T1*vcos(2*pi/365*t)
    tmps  = tmeds - T365*vcos(2*pi*(t - vint(t)))
    return (tmps, tmeds)
    
def Varianza(tmps, p):
#VARIANZA calcula la varianza de la distribucion del desarrollo 
#         en funcion de la temperatura
    numin = 0.000223
    nu = p[0] * tmps^2 + p[1] * tmps + p[2]
    #   nu = numin*ones(1,length(tmps))
    nu = (nu > numin) * nu + (nu <= numin) * numin 

def vonFoerster(dt, t, tau, nt, tmps, hmrs, pnu, fhnu, pdes, fhrates, pinput):
    """
    VONFOERSTER outputs evolution of an stage population based on vonFoster
    dt:      paso discreto de tiempo
	t:       vector de instantes de tiempos
	tau:     vector de instantes de tiempos
	nt:      numero total de espacios de tiempo
	tmps:    vector de temperaturas para todo instante discreto de tiempo
	hmrs:    vector de humedades relativas para todo instante discreto de tiempo
	pnu:     parametros de la funcion varianza()
	fhnu:    handler de la funcion varianza()
	pdes:    parametros de la funcion de desarrollo
	fhrates: handler de la funcion de desarrollo
	pinput:  vector de densidades de poblacion correspondiente al estado 
    """
    tol     = 0.0001                      # von Foerster model tolerance
    T,Tau   = np.meshgrid(t,tau)          # T,Tau matrices
    Pttau   = np.zeros((nt,nt))           # initialize matrix which will hold the convolution kernel
    wts     = np.zeros(nt)                # initialize vector which will hold weight for normalization
    rates   = np.zeros(nt)                # initialize vector of rates for each instant t
    hmrsnul = np.sum(abs(hmrs) > 0)
    #evalua con humedad y sin humedad
    if hmrsnul != 0:
       # calcula las tasas de desarrollo para temperaturas y humedades relativas dadas 
       for i in range(nt):
            rates[i] = feval(fhrates, tmps[i], hmrs[i], pdes)
    else:
        # calcula las tasas de desarrollo  para temperaturas dadas
        for i in range(nt):
            rates[i] = feval(fhrates, tmps[i], pdes) 
    
    RT = np.dot(dt * np.ones((nt,1)), [cumtrapz(rates)])  # create a matrix which is the cumulative development
    RTau  = np.transpose(RT)                                    # create transpose of cumulative matrix for use in kernel
    vexp  = np.vectorize(exp)   
    vsqrt = np.vectorize(sqrt)    
    
    if (pnu[0] == 0) and (pnu[1] == 0):
        nu = pnu[2]
        Pttau  = (T>Tau)*vexp(-(1-(RT-RTau))**2/(4*nu*(abs(T-Tau)+tol)))\
                           /vsqrt(4*pi*nu*(abs(T-Tau)**3+tol))         # extended von foerster kernel
    else:
        nus   = feval(fhnu, tmps, pnu)      # calcula las varianzas para temperatures dadas
        NU    = np.dot(np.ones(nt,1), nus)  # crea una matriz de varianzas en funcion de las temperaturas
        Pttau = (T>Tau)*exp(-(1-(RT-RTau))**2/(4*NU*(abs(T-Tau)+tol))) \
                           /vsqrt(4*pi*NU*(abs(T-Tau)**3+tol))         # extended von foerster kernel
    
    ints = dt*np.trapz(Pttau,axis=1)    # integrate in columns to normalize
    wts  = (ints>tol)*ints+(ints<=tol)  # calculate a  weighting factor 
                                        # make it one if the integral is too  small (< tol)
    pout = dt*np.dot(pinput/np.transpose(wts), Pttau)      # output distribution, normalized by integral of P
 
#Generacion de datos de prueba de una corrida 
startTime = timeit.default_timer()
tmin = 30         # dia inicial de la corrida 30 de enero
tmax = 30 + 100   # dia final de la corrida 30 de enero mas 100 dias
nd = tmax - tmin  # numero de dias
td = 4            # intervalos de tiempo diarios
nt = nd * td      # numero total de pasos de tiempo
dt = nd/float(nt)        # paso discreto de tiempo
t    = np.linspace(tmin,tmax-dt,nt) # vector de tiempos, t devuelve un vector de numeros entre un rango separados por el rango dividido la cantidad de intervalos
tau  = t                            # vector of tiempos, tau
T0   = 15 # T0   : temperatura media
T1   = 15 # T1   : amplitud termica anual
T365 = 15 # T365 : amplitud termica diaria 
tmps, tmeds = TempSim(t, T0, T1, T365)
hmrs  = np.zeros(nt)
# Datos para huevos de Diatrea 
fhrates = "BriereI"
pdes    = np.zeros(3) 
pdes[0] =  0.000131734
pdes[1] = 10.25740308  
pdes[2] = 36.65400490
# Descripcion de la funcion de varianza en funcion de las temperaturas
fhnu    = "Varianza"
pnu     = np.zeros(3) 
pnu[0]  = 0.0
pnu[1]  = 0.0  
pnu[2]  = 0.000223
# Datos de la poblacion de ingreso Pulso 100 individuos 
pinput  =  np.zeros(nt)
for i in range(4):
   pinput[i] = 25
idCorrida = "vFPy"
codfig = 8
vonFoerster(dt, t, tau, nt, tmps, hmrs, pnu, fhnu, pdes, fhrates, pinput)
endTime = timeit.default_timer()
print "Time elapsed: " + str((endTime - startTime)*1000) + " miliseconds."