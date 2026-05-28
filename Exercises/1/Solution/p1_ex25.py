# -*- coding: utf-8 -*-
"""
Created on Mon Feb  2 18:42:30 2026

@author: jfcte
"""
import numpy as np

# Python program to compute and print out the transmission and reflection probabilities using the formulas above.
# Data contained in the problem statement
# if we were to generalize this problem, all the program would only run if E>V, as if E<V all of the wave would end up being reflected

m = 9.11e-31 # mass of an electron (kg)
hbar = 1.0545718e-34 # in Js
eV = 1.602176634e-19 # 1 eV in J
E = 10*eV  # initial energy in J of the electron
V = 9*eV # height of the potential step in J

# Wavevectors
k1 = np.sqrt(2*m*E)/hbar #sqrt from numpy as np.sqrt
k2 = np.sqrt(2*m*(E-V))/hbar

# Probabilities T and R
T = (4*k1*k2)/(k1+k2)**2 # transmission
R = ((k1-k2)/(k1+k2))**2 # reflection

# Print the results
print(f"Transmision: T = {T:.4f}") # print the transmission probability with 4 decimals
print(f"Reflexion: R   = {R:.4f}") # print the reflection probability with 4 decimals
print(f"Comprobation: T + R = {T + R:.4f}") # print the sum of T+R, should be 1 if the code is well written
