# -*- coding: utf-8 -*-
from math import pi, sqrt
from numpy import zeros
from vpython import sphere, rate, color, canvas, vector

# Parámetros físicos
GMSol = 4 * pi**2

dt = 0.001
T = 76
tiempoTotal = T

N = int(tiempoTotal / dt)

# Vectores
x = zeros(N)
y = zeros(N)
r = zeros(N)
vx = zeros(N)
vy = zeros(N)
t = zeros(N)
v2 = zeros(N)

# Condiciones iniciales (perihelio)
t[0] = 0
x[0] = 0.59
y[0] = 0

vx[0] = 0
vy[0] = 11.3   # velocidad ajustada

v2[0] = vx[0]**2+vy[0]**2

r[0] = sqrt(x[0]**2 + y[0]**2)

# Integración Euler-Cromer
for i in range(N-1):

    vx[i+1] = vx[i] - GMSol * (x[i]/r[i]**3) * dt
    vy[i+1] = vy[i] - GMSol * (y[i]/r[i]**3) * dt

    x[i+1] = x[i] + vx[i+1] * dt
    y[i+1] = y[i] + vy[i+1] * dt

    r[i+1] = sqrt(x[i+1]**2 + y[i+1]**2)

    t[i+1] = t[i] + dt
    
    v2[i+1] = vx[i+1]**2+vy[i+1]**2


print("Distancia máxima (afelio):", max(r), "UA")
print("Velocidad máxima:", max(v2), "UA/año")


# Visualización
canvas(title="Órbita del Cometa Halley",
       width=700, height=700,
       center=vector(0,0,0),
       background=vector(0,0,0))

sol = sphere(pos=vector(0,0,0), color=color.yellow, radius=0.5)
cometa = sphere(pos=vector(x[0],y[0],0), color=color.white, radius=0.08)

for i in range(N-1):
    rate(200)

    px = x[i]
    py = y[i]

    sphere(pos=vector(px,py,0), radius=0.02, color=color.cyan)

    cometa.pos = vector(px,py,0)