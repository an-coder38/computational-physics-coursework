# -*- coding: utf-8 -*-
"""
Monte Carlo vectorizado (MUCHO más rápido)
"""

import numpy as np
from math import pi, gamma

# Parámetros
d = 12
R = 2
target_error = 5e-3 # para errores bajos hay error de RAM

# Volumen teórico
V_teorico = (pi**(d/2) / gamma(d/2 + 1)) * (R**d)
print("Volumen teórico:", V_teorico)

V_cubo = (2*R)**d

# Inicialización
NT = 100000
max_iter = 20

for k in range(max_iter):

    # Generar TODOS los puntos de golpe
    puntos = np.random.uniform(-R, R, size=(NT, d))

    # Suma de cuadrados por fila
    dist2 = np.sum(puntos**2, axis=1)

    # Contar puntos dentro
    NC = np.sum(dist2 <= R**2)

    # Estimación
    V_est = V_cubo * (NC / NT)
    error_rel = abs(V_est - V_teorico) / V_teorico

    print(f"N = {NT:>10} | V_est = {V_est:.6f} | error = {error_rel:.6e}")

    if error_rel < target_error:
        print("\n Precisión alcanzada")
        print("Número de números aleatorios necesarios ≈", NT * d)
        break

    NT *= 2

else:
    print("\n No se alcanzó la precisión deseada")
