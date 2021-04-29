# -*- coding: utf-8 -*-
"""
Física Computacional

Valentina Campos Aguilar
Luis Alfredo Guerrero Camacho

MÉTODO DIFERENCIAS FINITAS 
"""
#Se importan las bibliotecas necesarias para llevar a cabo el método escogido.
import numpy as np
import matplotlib.pyplot as plt

#Se definen los valores para las variables conocidas del problema a resolver.
D = 0.5
Lx = 10.0
Lt = 3.0
A = 2.0
x0 = 5.0
l = 1.5

''' Se crea la función que permite aproximar la densidad en los puntos (t,x) deseados en 
el material. 

Parámetros de entrada: matriz de densidades inicial, delta x, delta t, puntos (t,x) de la matriz.
Salida de la función: matriz de la densidad aproximada '''
def Aprox_Rho(rho_xti, dx, dt, puntos):
    #El orden de los ciclos indica que la densidad se obtiene iterando x a partir 
    #de la segunda columna para todo t. 
    for t in range(0, puntos[0]-1):
        for x in range(1, puntos[1]-1):
            #Se define la expresión matemática que permite obtener la densidad en el punto (t+1,x).
            rho_xti[t+1,x] = D*dt*(rho_xti[t,x+1]-2*rho_xti[t,x]+rho_xti[t,x-1])/dx**2 + rho_xti[t,x]
    return rho_xti


''' Se crea una función que permite obtener la solución de la aproximción de la densidad en
forma gráfica.

La misma no tiene parámetros de entrada.
Salida de la función: gráfico de contorno de la difusión, gráfico de suerficie de la densidad''' 
def Gráfica():
    #Se define el número de puntos de tiempo y de posición que contendrá la malla.
    tpuntos = 600
    xpuntos = 120
    #Se crean los arreglos de tiempo y posición con los límites y cantidad de puntos respectivos
    #de ambas magnitudes. 
    x = np.linspace(0, Lx, xpuntos)
    t = np.linspace(0, Lt, tpuntos)
    #Se crea la malla de dimensión posición x tiempo.
    X, T = np.meshgrid(x,t)
    #Se inicializa la matriz de dimensión(t,x) con los valores iniciales de la densidad. 
    rho_xti = np.zeros((tpuntos, xpuntos), float)
    #Se actualizan los valores de toda la primera fila de la matriz según la condición de 
    #contorno del problema.
    rho_xti[0,:] = A*np.exp(-(x-x0)**2/l)

    #Se definen delta x y delta t.
    dx = Lx/xpuntos
    dt = Lt/tpuntos
    
    #Se define la densidad aproximada llamando a la función Aprox_Rho con la matriz, los delta 
    #y delta t y los puntos(t,x) de la matriz anteriormente definidos, como parámetros de entrada. 
    P = Aprox_Rho(rho_xti, dx, dt, (tpuntos, xpuntos))
    #Se imprime la forma de la solución de la densidad para tener un control de las dimensiones.
    print(P.shape)
    
    #Se crea un gráfico de superficie que permite visualizar la aproximación de la densidad en 
    #función de la posición y del tiempo.
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.set_xlabel('x [m]')
    ax.set_ylabel('t [s]')
    ax.set_zlabel('ρ [densidad]')
    ax.plot_surface(X, T, P, rstride=1, cstride=1,
                cmap='cividis', edgecolor='none')
    ax.set_title('Gráfico de Superficie \n Aproximación de la densidad por Diferencias Finitas')
    plt.show()
    
    #Se crea un gráfico de contorno que permite visualizar los efectos de la posición y el 
    #tiempo en la difusión del material.
    plt.figure(figsize=(10,6))
    plt.scatter(X, T, c=P,cmap = 'jet')
    d = plt.colorbar()
    d.set_label('ρ [densidad]')
    plt.title('Gráfico de Contorno \n Comportamiento de la difusión por Diferencias Finitas')
    plt.xlabel('x [m]')
    plt.ylabel('t [s]')
    plt.show()
    
    
#Se llama a la función gráfica para la ejecución de todo el método y la obtención de los 
#resultados gráficos. 
Gráfica()

