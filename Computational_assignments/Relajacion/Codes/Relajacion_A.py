import numpy as np
import matplotlib.pyplot as plt
from numba import njit


# Función acelerada
@njit
def gauss_seidel(x, g, h, epsilon):
    N = len(x)
    delta = 1.0

    while delta > epsilon:
        delta = 0.0
        
        for i in range(1, N-1):
            old = x[i] # guardamos valor antiguo para obtener la diferencia entre ambos valores

            x[i] = 0.5 * (x[i-1] + x[i+1] + g*h**2) # fórmula sugerencia 2
            
            diff = abs(x[i] - old)
            if diff > delta:
                delta = diff

    return x, delta


# Parámetros
g = 9.81          # gravedad (m/s^2)
t0, tf = 0.0, 10.0
N = 500           # número de puntos
epsilon = 1e-6    # criterio de convergencia

# Malla temporal
t = np.linspace(t0, tf, N)
h = t[1] - t[0]

# Inicialización
x = np.zeros(N)

# Condiciones de contorno
x[0] = 0.0
x[-1] = 0.0

# Método de Gauss-Seidel
x, delta = gauss_seidel(x, g, h, epsilon)

print("Convergencia alcanzada:", delta)

# Solución analítica para comparar:
# x(t) = v0*t - 1/2 g t^2
# imponiendo x(10)=0 => v0 = g*10/2
v0 = g * tf / 2
x_exact = v0 * t - 0.5 * g * t**2

# Aproximación de la velocidad inicial numérica
v0_num = (x[1] - x[0]) / h

print(f"Velocidad inicial aproximada: {v0_num:.2f}")
print(f"Velocidad inicial exacta: {v0:.2f}")

# Gráfica
plt.figure()
plt.plot(t, x, label="Solución numérica")
plt.plot(t, x_exact, '--', label="Solución analítica")
plt.xlabel("Tiempo (s)")
plt.ylabel("Altura (m)")
plt.title("Trayectoria de la bola")
plt.legend()
plt.savefig('trayectoria_bola_A.png', dpi=150)
plt.show()

# Error numérico
error = np.abs(x - x_exact)

plt.figure()
plt.plot(t, error)
plt.xlabel("Tiempo (s)")
plt.ylabel("Error")
plt.title("Error numérico")
plt.savefig('errorA_2D.png', dpi=150)
plt.show()