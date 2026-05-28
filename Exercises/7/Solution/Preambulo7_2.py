import numpy as np
import matplotlib.pyplot as plt

# Parámetros
N = 7
J = 1.0
kB = 1.0
n_steps = 500

np.random.seed(38)

# Inicialización aleatoria
spins = np.random.choice([-1, 1], size=(N, N))

def energia(spins):
    E = 0
    for i in range(N):
        for j in range(N):
            s = spins[i, j]
            vecinos = (
                spins[(i+1)%N, j] +
                spins[i, (j+1)%N] +
                spins[(i-1)%N, j] +
                spins[i, (j-1)%N]
            )
            E += -J * s * vecinos * 0.5
    return E  # evitar doble conteo

def magnetizacion(spins):
    return np.sum(spins)

def metropolis(spins, T):
    for i in range(N):              # recorrido sistemático
        for j in range(N):
            s = spins[i, j]
            vecinos = (
                spins[(i+1)%N, j] +
                spins[i, (j+1)%N] +
                spins[(i-1)%N, j] +
                spins[i, (j-1)%N]
            )
            dE = 2 * J * s * vecinos

            if dE < 0 or np.random.rand() < np.exp(-dE/(kB*T)):
                spins[i, j] = -s

    return spins

# Barrido en temperatura
temperaturas = np.concatenate([
    np.linspace(1.0, 2.2, 10),     # antes de Tc
    np.linspace(2.2, 2.35, 10),    # alta resolución cerca de Tc
    np.linspace(2.35, 10, 30)     # después de Tc
])
E_T = []
M_T = []

for T in temperaturas:
    spins = np.random.choice([-1, 1], size=(N, N))

    # Equilibración
    for _ in range(1000):
        metropolis(spins, T)

    E_acum = 0
    M_acum = 0

    # Promedio
    for _ in range(n_steps):
        metropolis(spins, T)
        E_acum += energia(spins)
        M_acum += abs(magnetizacion(spins))

    E_T.append(E_acum / n_steps / (N*N))
    M_T.append(M_acum / n_steps / (N*N))

# Resultados
print("E(T):", E_T)
print("M(T):", M_T)

# Gráfica de energía
plt.figure()
plt.plot(temperaturas, E_T, 'o-')
plt.xlabel('T')
plt.ylabel('E(T)')
plt.title('Energía vs Temperatura')

# Gráfica de magnetización
plt.figure()
plt.plot(temperaturas, M_T, 'o-')
plt.xlabel('T')
plt.ylabel('M(T)')
plt.title('Magnetización vs Temperatura')
plt.show()