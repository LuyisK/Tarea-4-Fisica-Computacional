# -*- coding: utf-8 -*-
"""
Física Computacional

Valentina Campos Aguilar
Luis Alfredo Guerrero Camacho 

MÉTODO SERIES DE FOURIER
"""

#Se importan las bibliotecas necesarias para llevar a cabo el método escogido. 
import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt

#Se definen los valores para las variables conocidas del problema a resolver.
A = 2.0
D = 0.5
Lx = 10.0
Lt = 3.0
x0 = 5.0
l = 1.5

'''Se crea la siguiente función que permite obtener la densidad según la condición de contorno 
planteada en el problema.

Parámetro de entrada: x (posición)
Salida de la función: valor de la densidad descrito por la condición de contorno''' 
def Rho0(x):
    valorRho0 = A*np.exp(-(x-x0)**2/l)
    return valorRho0


''' Se define una función que permite calcular el coeficiente de Fourier. 

Parámetro de entrada: n(número de término)
Salida de la función: valor del coeficiente de Fourier '''
def Bn(n):
    #Se utiliza una función anónima para definir la función a integrar. 
    integral_bn = lambda x: (2/Lx)*Rho0(x)*np.sin(n*x*np.pi/Lx)
    #Se utiliza un método de integración de la biblioteca SciPy que permite 
    #calcular el valor del coeficiente con la función declarada anteriormente. 
    valorBn, error = integrate.quad(integral_bn, 0, Lx)
    return valorBn


'''Se crea la función que permite aproximar el valor de la densidad en el punto (x,t).

Parámetros de entrada: x(arreglo de posiciones), t (arreglo de tiempos), nt(número de términos)
Salida de la función: valor aproximado de la densidad'''
def AproxRho(x, t, nt):
    #Se inicializa el valor de la densidad en 0. 
    valorAproxRho = 0
    #Se crea un ciclo que realiza la sumatoria que corresponde al número de términos donde se 
    #aproxima la densidad para cada uno de estos. 
    for i in range(1, nt+1):
        #Se llama a la función Bn que toma el número de iteración como el parámetro de entrada
        #definiéndola como el valor del coeficiente. 
        coeficiente = Bn(i)
        #Se define la exresión matemática que permite hallar el valor aproximado de la densidad.
        valorAproxRho += coeficiente*np.sin(i*x*np.pi/Lx)*np.exp(-D*((i*np.pi/Lx)**2)*t)
    return valorAproxRho


#Se define el número de puntos de tiempo y de posición que contendrá la malla. 
xpuntos = 120
tpuntos = 600

#Se crean los arreglos de posición y tiempo con los límites y cantidad de puntos respectivos
#de ambas magnitudes. 
x = np.linspace(0, Lx, xpuntos) 
t = np.linspace(0, Lt, tpuntos)

#Se crea la malla de dimensión posición x tiempo. 
X, T = np.meshgrid(x,t)

#Se define el número de términos a calcular. 
númerotérminos = 5

#Se define el la densidad aproximada llamando a la función AproxRho con el número de
#términos y la malla establecidos como parámetros.
P = AproxRho(X, T, númerotérminos)
#Se imprime la forma de la solución de la densidad para tener un control de las dimensiones. 
print(P.shape)

'''Se crea un gráfico de superficie que permite visualizar la aproximación de la densidad en 
función de la posición y del tiempo'''
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.set_xlabel('x [m]')
ax.set_ylabel('t [s]')
ax.set_zlabel('ρ [densidad]')
ax.plot_surface(X, T, P, rstride=1, cstride=1,
                cmap='cividis', edgecolor='none')
ax.set_title('Gráfico de Superficie \n Aproximación de la densidad por Series de Fourier')
plt.show()

'''Se crea un gráfico de contorno que permite visualizar los efectos de la posición y el 
tiempo en la difusión del material'''
plt.figure(figsize=(10,6))
plt.scatter(X, T, c=P,cmap = 'jet')
d = plt.colorbar()
d.set_label('ρ [densidad]')
plt.title('Gráfico de Contorno \n Comportamiento de la difusión por Series de Fourier')
plt.xlabel('x [m]')
plt.ylabel('t [s]')
plt.show()

    
    