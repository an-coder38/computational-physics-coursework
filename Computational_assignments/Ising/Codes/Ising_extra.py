# Código que estudia de forma más detallada el modelo de Ising 2D
import numpy as np
import matplotlib.pyplot as plt
from numba import njit

# Parámetros

N = 8       # tamaño de la red (N x N)
J = 1.0     # constante de intercambio ferromagnético
kB = 1.0    # constante de Boltzmann (unidades reducidas)

temperaturas = [0.25, 1, 2, 3.5, 10, 50]

H_max = 4
dH = 0.25

# Barridos de campo: subida de -H_max a +H_max y bajada de vuelta
H_values_up   = np.arange(-H_max, H_max + dH, dH)
H_values_down = np.arange( H_max, -H_max - dH, -dH)

# n_eq:   pasos de Metropolis entre medidas (equilibración local)
# n_meas: número de medidas
n_eq   = 50
n_meas = 2000

np.random.seed(35)  # semilla global para numpy fuera de numba

# Funciones utilizadas y aceleradas con numba

@njit
def fijar_seed(seed):
    """Fija la semilla de numpy dentro del entorno numba para reproducibilidad."""
    np.random.seed(seed)

@njit
def magnetizacion(spins):
    """Magnetización media por spin: M = sum(spins) / (N*N)"""
    return np.sum(spins) / (N * N)

@njit
def energia(spins, H):
    """
    Energía total del sistema Ising con campo externo H.
    Hamiltoniano: E = -J * sum(s_i * s_j) - H * sum(s_i)
    El factor 1/2 evita contar cada par de vecinos dos veces.
    """
    E = 0.0
    for i in range(N):
        for j in range(N):
            s = spins[i, j]
            # Vecinos con condiciones de contorno periódicas
            vecinos = (
                spins[(i+1) % N, j] +
                spins[i, (j+1) % N] +
                spins[(i-1) % N, j] +
                spins[i, (j-1) % N]
            )
            E += -J * s * vecinos - H * s
    return E / 2.0

@njit
def metropolis(spins, H, T):
    """
    Un paso de barrido Metropolis sobre toda la red.
    Para cada spin calcula dE al flippear y acepta si:
      - dE < 0 (energéticamente favorable), o
      - con probabilidad exp(-dE / kB*T) (fluctuación térmica)
    """
    for i in range(N):
        for j in range(N):
            s = spins[i, j]
            vecinos = (
                spins[(i+1) % N, j] +
                spins[i, (j+1) % N] +
                spins[(i-1) % N, j] +
                spins[i, (j-1) % N]
            )
            dE = 2.0 * s * (J * vecinos + H)  # variación de energía al flipear s_i

            if dE < 0.0 or np.random.random() < np.exp(-dE / (kB * T)):
                spins[i, j] = -s

# Simulación H

fijar_seed(35)

resultados = {}

for T in temperaturas:
    print(f"Simulando T = {T}...")

    # Estado inicial: todos los spins alineados (ferromagnético)
    spins = np.ones((N, N), dtype=np.int32)

    M_up, M_down = [], []
    E_up, E_down = [], []

    # --- SUBIDA: H de -H_max a +H_max ---
    for H in H_values_up:
        m_sum = e_sum = 0.0

        for _ in range(n_meas): # n_meas medidas
            # n_eq pasos de Metropolis entre medidas para decorrelacionar
            for _ in range(n_eq):
                metropolis(spins, H, T)

            m_sum += magnetizacion(spins)
            e_sum += energia(spins, H)

        M_up.append(m_sum / n_meas) # tomamos la media de las medidas
        E_up.append(e_sum / n_meas)

    # --- BAJADA: H de +H_max a -H_max ---
    # Se parte del estado final de la subida para ver si hay histéresis
    for H in H_values_down:
        # Equilibración inicial al nuevo H antes de medir
        for _ in range(n_eq):
            metropolis(spins, H, T)

        m_sum = e_sum = 0.0

        for _ in range(n_meas):
            metropolis(spins, H, T)
            m_sum += magnetizacion(spins)
            e_sum += energia(spins, H)

        M_down.append(m_sum / n_meas)
        E_down.append(e_sum / n_meas)

    # guardar los resultados 
    resultados[T] = {
        "H_up": H_values_up, "M_up": M_up, "E_up": E_up,
        "H_down": H_values_down, "M_down": M_down, "E_down": E_down
    }


# Barrido de temperatura
H_fijos  = np.arange(-4,6,2)
T_values = np.arange(0.4, 10.2, 0.2)

# Simulación E(T) para H fijo

resultados_T = {}

for H in H_fijos:
    print(f"Simulando E(T) con H = {H}...")
    E_vs_T = []

    # Estado inicial ferromagnético; se arrastra entre temperaturas
    spins = np.ones((N, N), dtype=np.int32)

    for T in T_values:
        e_sum = 0.0

        for _ in range(n_meas):
            for _ in range(n_eq):
                metropolis(spins, H, T)
            e_sum += energia(spins, H)

        E_vs_T.append(e_sum / n_meas)

    resultados_T[H] = np.array(E_vs_T)

# Gráficas M

# --- Global (todas las temperaturas en una gráfica) ---
plt.figure(figsize=(8, 6))
for T in temperaturas:
    data = resultados[T]
    H = np.concatenate([data["H_up"], data["H_down"]]) # unimos H_up y -h_down para luego graficarlas
    M = np.concatenate([data["M_up"], data["M_down"]])
    plt.plot(H, M, label=f"T={T}")
plt.xlabel("H")
plt.ylabel("M")
plt.legend()
plt.grid()
plt.savefig("magnetizacion_hist_final.png", dpi=300, bbox_inches="tight")
plt.show()

# --- Subplots por temperatura (una gráfica por temperatura) pero en una misma imagen ---
fig, axes = plt.subplots(1, len(temperaturas), figsize=(18, 5), sharey=True)
for idx_T, T in enumerate(temperaturas):
    ax = axes[idx_T]
    data = resultados[T]
    ax.plot(data["H_up"],   data["M_up"],   '-',  label='Subida')
    ax.plot(data["H_down"], data["M_down"], '--', label='Bajada')
    ax.set_title(f'T = {T}')
    ax.set_xlabel('H')
    ax.grid()
    ax.legend()
    if idx_T == 0:
        ax.set_ylabel('Magnetización M')
plt.tight_layout()
plt.savefig("histeresis_vs_temperatura_subplots.png", dpi=300, bbox_inches='tight')
plt.show()

# Gráficas E(H) para T fijo

# --- Global (subida y bajada diferenciadas) ---
plt.figure(figsize=(8, 6))
for T in temperaturas:
    data = resultados[T]
    plt.plot(data["H_up"],   data["E_up"],   '-',  label=f"T={T} subida")
    plt.plot(data["H_down"], data["E_down"], '--', label=f"T={T} bajada")
plt.xlabel("H")
plt.ylabel("Energía")
plt.legend()
plt.savefig("energia_hist_final.png", dpi=300, bbox_inches="tight")
plt.show()

# --- Subplots por temperatura ---
fig, axes = plt.subplots(1, len(temperaturas), figsize=(18, 5), sharey=True)
for idx_T, T in enumerate(temperaturas):
    ax = axes[idx_T]
    data = resultados[T]
    ax.plot(data["H_up"],   data["E_up"],   '-',  label='Subida')
    ax.plot(data["H_down"], data["E_down"], '--', label='Bajada')
    ax.set_title(f'T = {T}')
    ax.set_xlabel('H')
    ax.grid()
    ax.legend()
    if idx_T == 0:
        ax.set_ylabel('Energía E')
plt.tight_layout()
plt.savefig("energia_vs_temperatura_subplots.png", dpi=300, bbox_inches='tight')
plt.show()

# Gráficas E(T) para H fijo

# --- Global ---
plt.figure(figsize=(8, 6))
for H in H_fijos:
    plt.plot(T_values, resultados_T[H], label=f"H={H}")
plt.xlabel("Temperatura T")
plt.ylabel("Energía E")
plt.title("Energía vs Temperatura (H fijo)")
plt.legend()
plt.grid()
plt.savefig("energia_vs_T_global_zoom.png", dpi=300, bbox_inches="tight")
plt.show()

# --- Subplots por H ---
fig, axes = plt.subplots(1, len(H_fijos), figsize=(15, 5), sharey=True)
for idx_H, H in enumerate(H_fijos):
    ax = axes[idx_H]
    ax.plot(T_values, resultados_T[H], 'o-')
    ax.set_title(f"H = {H}")
    ax.set_xlabel("T")
    ax.grid()
    if idx_H == 0:
        ax.set_ylabel("Energía E")
plt.suptitle("Energía vs Temperatura por H")
plt.tight_layout()
plt.savefig("energia_vs_T_subplots_zoom.png", dpi=300, bbox_inches="tight")
plt.show()