import numpy as np
import matplotlib.pyplot as plt

# Parámetros
D = 1
L = 50

dx = 1
dy = 1

dt = 1 / (2*D*(1/dx**2 + 1/dy**2))  # estabilidad para dx no necesariamente igual a dy

NX = int(2*L/dx) + 1
NY = int(2*L/dy) + 1
NT = 200

# Mallas
x = np.linspace(-L, L, NX)
y = np.linspace(-L, L, NY)

# Densidad
rho = np.zeros((NX, NY, NT+1))

# # Condición inicial: pico en el centro
# cx, cy = NX//2, NY//2
# rho[cx, cy, 0] = 1

# Condición inicial para evitar que se distinga la paridad de los puntos: bloque [-5,5]x[-5,5]
for i in range(NX):
    for j in range(NY):
        if -5 <= x[i] <= 5 and -5 <= y[j] <= 5:
            rho[i, j, 0] = 1

# Evolución temporal sin vectorizar
# for n in range(NT):
#     for i in range(1, NX-1):
#         for j in range(1, NY-1):
#             rho[i, j, n+1] = rho[i, j, n] + D*dt*(
#                 (rho[i+1, j, n] - 2*rho[i, j, n] + rho[i-1, j, n]) / dx**2 +
#                 (rho[i, j+1, n] - 2*rho[i, j, n] + rho[i, j-1, n]) / dy**2
#             )

# Evolución temporal vectorizada
for n in range(NT):
    rho[1:-1,1:-1,n+1] = rho[1:-1,1:-1,n] + D*dt*(
        (rho[2:,1:-1,n] - 2*rho[1:-1,1:-1,n] + rho[:-2,1:-1,n]) / dx**2 +
        (rho[1:-1,2:,n] - 2*rho[1:-1,1:-1,n] + rho[1:-1,:-2,n]) / dy**2
    )

# --- REPRESENTACIÓN ---

def plot_estado(paso, titulo):
    plt.figure()
    plt.imshow(rho[:, :, paso], extent=[-L, L, -L, L], origin='lower', aspect='auto', cmap='jet')
    plt.colorbar(label=r"$\rho$")
    plt.title(titulo)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.show()

# Instantes pedidos
plot_estado(0, "t = 0")
plot_estado(50, "t = 50")
plot_estado(200, "t = 200")