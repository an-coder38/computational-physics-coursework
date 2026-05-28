import numpy as np
import matplotlib.pyplot as plt

# Parámetros
m = 2500     # caminantes
N = 2000    # pasos

np.random.seed(10) # por reproducibilidad

# bloque que comprueba cómo funcionan las condiciones periódicas de frontera para L=35
L = 50
for val in [-35, 35, 51, -51, 120, -120]:
    x = (val + L) % (2*L) - L
    print(val, "pasa a", x)

# Inicia con distribución uniforme en [-5,5]x[-5,5]
x = np.random.uniform(-5, 5, size=m)
y = np.random.uniform(-5, 5, size=m)

# Entropía
S = []

# Malla fija
L = 50
bins = 10
edges = np.linspace(-L, L, bins+1) #creo dimensiones para el futuro histograma con la división por celdas

for n in range(N):

    # Movimiento aleatorio
    pasos_x = np.random.choice([-1, 1], size=m)
    pasos_y = np.random.choice([-1, 1], size=m)

    x += pasos_x
    y += pasos_y
    
    x = (x + L) % (2*L) - L # condiciones periódicas en x
    y = (y + L) % (2*L) - L # condiciones periódicas en y

    # Histograma 2D dividido en las celdas definidas anteriormente
    H, _, _ = np.histogram2d(x, y, bins=[edges, edges])

    # renormalización
    if np.sum(H) > 0:
        P = H / np.sum(H)
    else:
        P = H

    # Evitar log(0)
    P_nonzero = P[P > 0]

    # Entropía
    S_t = -np.sum(P_nonzero * np.log(P_nonzero))
    S.append(S_t)

S = np.array(S)
t = np.arange(N)

# Gráfica de la entropía

plt.figure(figsize=(8,5))
plt.plot(t, S, label="S(t) paseo aleatorio")

# Entropía máxima teórica
S_max = np.log(bins * bins)
plt.axhline(S_max, linestyle='--', label="S máxima teórica = ln(100)")
plt.xlabel("tiempo / pasos")
plt.ylabel("S(t)")
# plt.title("Evolución de la entropía")
plt.legend()
plt.grid()
plt.savefig("evolucion_entropia.png", dpi=300)
plt.show()
