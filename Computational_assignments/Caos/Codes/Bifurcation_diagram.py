import numpy as np
import matplotlib.pyplot as plt

def pendulo_no_lineal(theta0, omega0, timeTotal, t0, dt, q, FI, OmegaI):

    N = int((timeTotal - t0) / dt)
    theta = np.zeros(N)
    omega = np.zeros(N)
    t = np.zeros(N)
    theta[0] = theta0
    omega[0] = omega0
    t[0] = t0

    for i in range(N-1):
        # Aceleración angular (usando g/l = 1)
        a = -np.sin(theta[i]) - q * omega[i] + FI * np.sin(OmegaI * t[i])
        omega[i+1] = omega[i] + a * dt
        theta[i+1] = theta[i] + omega[i+1] * dt
        t[i+1] = t[i] + dt

    return t, theta

# Condiciones iniciales, las mismas que en P4
g = 9.8
l = 9.8
theta0 = 0.20
omega0 = 0.0
t0 = 0.0
dt = 0.01

q = 0.5
OmegaI = 2/3

# Período de la fuerza externa, en sus múltiplos tomaremos cuenta de la posición angular del sistema
T = 2 * np.pi / OmegaI

# Se realiza un barrido cerca del valor de FI (se realizó previamente algún barrido más extenso)
FI_min = 1.1
FI_max = 1.6 
num_FI = 100 # número de valores de FI
valores_FI = np.linspace(FI_min, FI_max, num_FI)

n_trans = 20 # número de períodos para el transitorio
n_puntos = 50 # número de puntos a guardar por cada FI

timeTotal = (n_trans + n_puntos) * T # tiempo total de simulación

FI_plot = [] # listas inicializadas
theta_plot = []

# Bucle sobre los valores de FI
print("Diagrama de bifurcación")
for FI in valores_FI:
    # print('FI')
    t, theta = pendulo_no_lineal(theta0, omega0, timeTotal, t0, dt, q, FI, OmegaI)

    # Determinar los índices correspondientes a los tiempos de Poincaré
    # Después del transitorio, queremos t = n * T para n = n_trans+1 hasta n_trans+n_puntos
    indices_poincare = []
    for n in range(n_trans + 1, n_trans + n_puntos + 1):
        t_poinc = n * T
        # Encontrar el índice más cercano en el array de tiempos
        idx = int(round((t_poinc - t0) / dt))
        # Verificar que esté dentro del rango (por seguridad)
        if 0 <= idx < len(t):
            indices_poincare.append(idx)

    # Extraer los ángulos en esos índices 
    for idx in indices_poincare:
        theta_mod = theta[idx] % (2*np.pi) # aplicar theta < 2*pi
        FI_plot.append(FI)
        theta_plot.append(theta_mod)

# Figura del diagrama de bifurcación
plt.figure(figsize=(10, 6))
plt.scatter(FI_plot, theta_plot, c='k', marker='.')
plt.xlabel(r'F_I')
plt.ylabel(r'$\theta$ (rad)')
plt.title('Diagrama de bifurcación del péndulo forzado:/n $q=0.5,\ \Omega=2/3,\ \theta_0=0.20,\ \omega_0=0$')
plt.xlim(FI_min, FI_max)
plt.ylim(0, 2*np.pi)
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("bif.png", dpi=300)
plt.show()

#zoom al diagrama de bifurcación en el que se ve el doblado del periodo y luego el caos
plt.figure(figsize=(10, 6))
plt.scatter(FI_plot, theta_plot, c='k', marker='.')
plt.xlabel(r'F_I')
plt.ylabel(r'$\theta$ (rad)')
plt.title('Diagrama de bifurcación del péndulo forzado:/n $q=0.5,\ \Omega=2/3,\ \theta_0=0.20,\ \omega_0=0$')
plt.xlim(1.32, 1.50)
plt.ylim(0, 4)
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("bif_zoom.png", dpi=300)
plt.show()
