# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

# Parámetros
N = 8
J = 1.0
kB = 1.0
T = 0.25

np.random.seed(38)

# Inicialización: todos alineados (importante para histéresis)
spins = np.ones((N, N), dtype=int)

# Función de magnetización
def magnetizacion(spins):
    return np.sum(spins) / (N*N)

# Paso de Metropolis con campo H
def metropolis(spins, H):
    for i in range(N):                 # recorrido sistemático
        for j in range(N):
            s = spins[i, j]
            vecinos = (
                spins[(i+1)%N, j] +
                spins[i, (j+1)%N] +
                spins[(i-1)%N, j] +
                spins[i, (j-1)%N]
            )
            dE = 2 * s * (J * vecinos + H)

            if dE < 0 or np.random.rand() < np.exp(-dE/(kB*T)):
                spins[i, j] = -s

    return spins

# Barrido de campo (ida y vuelta)
H_max = 4
dH = 0.5

H_values_up = np.arange(-H_max, H_max + dH, dH)
H_values_down = np.arange(H_max, -H_max - dH, -dH)

# Diferentes tiempos Monte Carlo
mc_steps_list = [10, 100, 1000]

resultados = {}

for mc_steps in mc_steps_list:
    spins = np.ones((N, N), dtype=int)

    M_up = []
    M_down = []

    # Barrido ascendente
    for H in H_values_up:
        for _ in range(mc_steps):
            metropolis(spins, H)
        M_up.append(magnetizacion(spins))

    # Barrido descendente
    for H in H_values_down:
        for _ in range(mc_steps):
            metropolis(spins, H)
        M_down.append(magnetizacion(spins))

    H_total = np.concatenate((H_values_up, H_values_down))
    M_total = np.concatenate((M_up, M_down))

    # Guardar resultados sin necesidad de crear un nuevo archivo .txt o .csv
    resultados[mc_steps] = {
        "H_up": H_values_up,
        "M_up": M_up,
        "H_down": H_values_down,
        "M_down": M_down,
        "H_total": H_total,
        "M_total": M_total
    }

# Gráfica conjunta

plt.figure()

for mc_steps in mc_steps_list:
    data = resultados[mc_steps]

    H_total = np.concatenate((data["H_up"], data["H_down"]))
    M_total = np.concatenate((data["M_up"], data["M_down"]))

    plt.plot(H_total, M_total, 'o-', label=f'MC steps = {mc_steps}')

plt.xlabel('Campo externo H')
plt.ylabel('Magnetización M')
plt.title('Histéresis (curvas superpuestas)')
plt.legend()
plt.grid()

# Subplots

fig, axes = plt.subplots(1, len(mc_steps_list), figsize=(15, 4), sharey=True)

for idx, mc_steps in enumerate(mc_steps_list):
    data = resultados[mc_steps]
    ax = axes[idx]

    ax.plot(data["H_up"], data["M_up"], 'o-', label='Subida')
    ax.plot(data["H_down"], data["M_down"], 's--', label='Bajada')

    ax.set_title(f'MC steps = {mc_steps}')
    ax.set_xlabel('H')
    ax.grid()
    ax.legend()

    if idx == 0:
        ax.set_ylabel('Magnetización M')

plt.suptitle('Histéresis en el modelo de Ising (8x8, T=0.25)')
plt.tight_layout()

plt.show()