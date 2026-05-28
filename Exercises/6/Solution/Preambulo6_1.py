# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 14:41:15 2026

@author: jfcte
"""

# Llamada a funciones externas
import numpy as np
from random import random, seed
from pylab import plot, show, figure, xlim, xlabel, ylabel, title

zeros = np.zeros

# Inicialización de variables, tablas y listas
n = 100 # col, número de pasos
k=200 # filas, número de caminantes
x = zeros((n,k))
x2 = zeros((n,k))
x[0,0:k-1] = 0
x2[0,0:k-1] = x[0,0:k-1]**2
paso = 0

seed(10) # Semilla del generador

# Bucle principal para los "n" pasos
for j in range (k):    
    for i in range(n-1):
        r = random()            # Se escoge número aleatorio
        if r < 0.5:             # Condición
            paso = 1            # Suma un paso hacia +x
        else:
            paso = -1           # Suma un paso hacia -x
        x[i+1,j] = x[i,j] + paso    # Posición después de i+1 pasos
        x2[i+1,j] = x[i+1,j]**2     # Cuadrado de la posición anterior
        # print(r, " ", paso, " ", x[i+1], " ", x2[i+1])
        
# Representación de las posiciones al cuadrado en función de los pasos
promedio = np.mean(x2, axis=1)
figure(2)
plot(promedio, 'o')
xlim(0)
xlabel("Pasos")
ylabel(r"$x^{2}$")
title('Paseo aleatorio de un caminante en una dimensión')
show()

# Representación de las posiciones en función de los pasos
figure(1)    
for j in range (k):
    plot(x[:,j], "o")
xlabel("Pasos")
ylabel(r"$x$")
title('Paseo aleatorio de un caminante en una dimensión')
xlim(0)
show()