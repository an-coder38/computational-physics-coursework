# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 14:59:28 2026

@author: jfcte
"""

import numpy as np
import matplotlib.pyplot as plt

# Paseo aleatorio

m = 1000   # caminantes
N = 100    # pasos

np.random.seed(8)

# Inicia con distribución uniforme en [-5,5]x[-5,5]
x = np.random.uniform(-5, 5, size=m)
y = np.random.uniform(-5, 5, size=m)

msd_rw = []

for n in range(N):
    pasos_x = np.random.choice([-1, 1], size=m) # escoge +1 o -1 de forma aleatoria con misma probabilidad, más rápido que otro bucle for
    pasos_y = np.random.choice([-1, 1], size=m)

    x += pasos_x
    y += pasos_y

    r2 = x**2 + y**2
    msd_rw.append(np.mean(r2))

msd_rw = np.array(msd_rw)
t_rw = np.arange(N)

# Difusión en 2D

D = 0.5   # valor esperado teórico
L = 50
dx = dy = 1

dt = 1 / (2 * D * (1/dx**2 + 1/dy**2))  # condición de estabilidad que se cumple si dx y dy no son iguales

NX = int(2*L/dx) + 1
NY = int(2*L/dy) + 1
NT = 200

x_grid = np.linspace(-L, L, NX)
y_grid = np.linspace(-L, L, NY)

rho = np.zeros((NX, NY, NT+1))

# Condición inicial (bloque centrado en 0,0 de tamaño 5^2)
for i in range(NX):
    for j in range(NY):
        if -5 <= x_grid[i] <= 5 and -5 <= y_grid[j] <= 5:
            rho[i, j, 0] = 1

# Evolución temporal (vectorizada)
for n in range(NT):
    rho[1:-1,1:-1,n+1] = rho[1:-1,1:-1,n] + D*dt*(
        (rho[2:,1:-1,n] - 2*rho[1:-1,1:-1,n] + rho[:-2,1:-1,n]) / dx**2 +
        (rho[1:-1,2:,n] - 2*rho[1:-1,1:-1,n] + rho[1:-1,:-2,n]) / dy**2)

#crea matrices de X e Y para poder utilizar la distribución de probabilidad en xi, yj
X, Y = np.meshgrid(x_grid, y_grid, indexing='ij') # x en la que cada fila tiene el mismo valor, y tiene en cada columna el mismo valor, así puedo localizar en ij las posiciones y luego aplicar la distribución de probabilidad

#necesario para <r^2> = sum_{i,j} (x_i^2 + y_j^2) * rho_ij / sum_{i,j} rho_ij
msd_diff = [] #inicializo la lista en la que guardar los resultados

for n in range(NT+1):
    r2 = X**2 + Y**2
    media = np.sum(r2 * rho[:, :, n]) / np.sum(rho[:, :, n]) # cálculo del valor esperado de la distancia al cuadrado
    msd_diff.append(media)

msd_diff = np.array(msd_diff)
t_diff = np.arange(NT+1) * dt

# Ajusto los datos a una recta
coef_rw = np.polyfit(t_rw, msd_rw, 1)
coef_diff = np.polyfit(t_diff, msd_diff, 1)

print("=== RESULTADOS ===")
print(f"Paseo aleatorio: pendiente ≈ {coef_rw[0]:.3f} (esperado ~2.0)")
print(f"Difusión: pendiente ≈ {coef_diff[0]:.3f} (esperado ~4D = {4*D:.1f})")

# Plots de comparación rw con difusión
plt.figure(figsize=(8,5))
plt.plot(t_rw, msd_rw, label="Random Walk 2D") # paseo aleatorio
plt.plot(t_diff, msd_diff, label="Difusión 2D") #difusión
plt.xlabel("tiempo / pasos")
plt.ylabel("<r²>")
plt.title("Comparación: Paseo aleatorio vs Difusión")
plt.legend()
plt.grid()
plt.savefig("rw_vs_dif.png", dpi=300)
plt.show()

# --- REPRESENTACIÓN ---

def plot_estado(paso, titulo):
    # representa la distribución de probabilidad para diferentes instantes temporales
    plt.figure()
    plt.imshow(rho[:, :, paso], extent=[-L, L, -L, L], origin='lower', aspect='auto', cmap='jet')
    plt.colorbar(label=r"$\rho$")
    plt.title(titulo)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.savefig(f"t={paso}.png", dpi=300)
    plt.show()

plot_estado(0, "t = 0")
plot_estado(NT//4, f"t = {NT//4}")
plot_estado(NT, f"t = {NT}")
