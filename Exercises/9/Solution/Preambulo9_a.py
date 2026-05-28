import numpy as np
import matplotlib.pyplot as plt
from numba import njit

# Función acelerada
@njit
def poisson_solver(phi, rho, epsilon):
    N = phi.shape[0]
    delta = 1.0
    
    while delta > epsilon:
        delta = 0.0
        
        for i in range(1, N-1):
            for j in range(1, N-1):
                old_phi = phi[i, j]
                
                phi[i, j] = 0.25 * (
                    phi[i+1, j] + phi[i-1, j] +
                    phi[i, j+1] + phi[i, j-1] +
                    rho[i, j]
                )
                
                diff = abs(phi[i, j] - old_phi)
                if diff > delta:
                    delta = diff

    return phi, delta

# Parámetros
N = 100
epsilon = 1e-6
phi = np.zeros((N, N))
rho = np.zeros((N, N))

rho[60:80, 60:80] = +1.0
rho[20:40, 20:40] = -1.0

# Ejecutar
phi, delta = poisson_solver(phi, rho, epsilon)

print("Convergencia alcanzada:", delta)

plt.imshow(phi, origin='lower')
plt.colorbar(label="Potencial")
plt.title("Solución de la ecuación de Poisson")
plt.show()