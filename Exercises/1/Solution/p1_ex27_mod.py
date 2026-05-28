# -*- coding: utf-8 -*-
"""
Created on Mon Feb  2 18:42:30 2026

@author: jfcte
"""
# Write a program that prints in increasing order all Catalan numbers
import numpy as np
import pylab as pyl
from matplotlib.ticker import MaxNLocator

print("Este es un programa que calcula los números catalanes")
# It keeps the limit as a float because the limit is not necessary an integer (as Catalan numbers are)
lim=float(input("Escriba el límite superior (no incluido) hasta el que calcularlos:"))
c=1 #first catalan number
n=0
c_exact=[] #c obtained with the recurrence sequence
c_estim=[] #c aproximated with c_rec=4^n/((n^1.5)*sqrt(pi))
while c<lim: #notice that the last calculated catalan number may be calculated but not printed
    c_exact.append(c)
    if n>0:
        c_estim.append((4**n)/((n**1.5)*np.sqrt(np.pi)))
    else:
        c_estim.append(0.0) # 0 so it avoids a divergence
    print(int(c)) #prints c as an integer as the catalan numbers are integers
    c=(4*n+2)*c/(n+2) #defining c(n+1) with c(n)
    n+=1
print(f"Hay {n} números catalanes inferiores a {lim}")

# a comparison between the catalan numbers obtained via the recurrence sequence and its aproximation can be graphed

n_vals = range(1, len(c_exact) + 1) # the list of the catalan numbers, it starts in 1 instead of 0 so the fist catalan number will be paired and graphed with n=1

# Gráfica
pyl.plot(n_vals, c_exact, "bo-",
        label="Recurrencia",
        markeredgecolor="black",
        markeredgewidth=0.5)

pyl.plot(n_vals, c_estim, "rs--",
        label="Estimación asintótica",
        markeredgecolor="black",
        markeredgewidth=0.5)

pyl.gca().xaxis.set_major_locator(MaxNLocator(integer=True)) # so the x axis only shows integer numbers
# from pylab, get current axes, select x_axis, select the principal units of the axis, only integers as principal marks in the axis

pyl.xlabel("n")
pyl.ylabel("Catalan number")
pyl.title("Números Catalanes: exactos vs estimación")

pyl.legend()
pyl.grid(True)
pyl.show()
