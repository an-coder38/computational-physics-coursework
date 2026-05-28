# condensador con cilindro interno a 2.5V y externo a tierra, las tapas en 0 y L están a un potencial cualquiera -- problema independiente de z y theta
import numpy as np
import matplotlib.pyplot as plt
from numba import njit

# Parámetros
Ri, Re, L = 0.01, 0.04, 0.20 # radios interno, externo y longitud (L no se usa en 2D)
Vi, Ve    = 2.5, 0.0 # potencial en el cilindro interno y externo

# Malla radial
N = 200 # número de puntos en r
r = np.linspace(Ri, Re, N) # malla radial uniforme
h = r[1] - r[0] # paso radial

# Inicializamos el potencial
V = np.linspace(Vi, Ve, N) # condición inicial lineal ("sentido común")
epsilon = 1e-10 # criterio de convergencia
max_iter = 100000 # máximo de iteraciones

# Gauss-Seidel en cilíndricas
@njit
def gauss_seidel_cilindrico(V, r, h, Vi, Ve, epsilon, max_iter):

    N = len(V) # número de nodos
    errores = np.zeros(max_iter) # para guardar el error

    for it in range(max_iter): # evitar colapsos

        delta = 0.0 # error máximo de esta iteración

        for i in range(1, N - 1): # puntos interiores

            old = V[i] # guardamos valor anterior para el error

            # factor geométrico cilíndrico
            fac = h / (2.0 * r[i])

            # actualización Gauss-Seidel
            V[i] = (V[i + 1] * (1.0 + fac) + V[i - 1] * (1.0 - fac)) / 2.0

            # error
            diff = abs(V[i] - old)

            # máximo error global
            if diff > delta:
                delta = diff

        # condiciones de contorno
        V[0] = Vi
        V[-1] = Ve

        # guardamos convergencia
        errores[it] = delta

        # criterio de break
        if delta < epsilon:
            return V, errores[:it+1], it+1

    return V, errores, max_iter


# Corremos la función que resuelve el problema
V, errores, it_conv = gauss_seidel_cilindrico(V, r, h, Vi, Ve, epsilon, max_iter)

print(f"Convergencia en {it_conv} iteraciones (error={errores[-1]:.2e})")

# Solución analítica para comparar
V_an = Vi * np.log(Re / r) / np.log(Re / Ri)

# # Tabla que compara analítico con numérico
# print(f"\n Ri={Ri*100} cm, Re={Re*100} cm, L={L*100} cm, Vi={Vi} V, Ve={Ve} V")
# print(f" {'r [cm]':>8} {'V_GS':>10} {'V_anal':>10} {'|err|':>10}")

# for i in [0, N//4, N//2, 3*N//4, N-1]:
#     print(f"{r[i]*100:8.3f} {V[i]:10.6f} {V_anal[i]:10.6f} {abs(V[i]-V_anal[i]):10.2e}")

# Solución radial V(r)
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

# V(r)
ax = axes[0]
ax.plot(r*100, V, lw=2.5, label='Gauss-Seidel (numérico)')
ax.plot(r*100, V_an, lw=2, ls='--', label=r'Analítico')
ax.set_xlabel('r [cm]')
ax.set_ylabel('V(r) [V]')
ax.set_title('Potencial en el cilindro')
ax.legend()
ax.set_xlim(Ri*100, Re*100)

# grafico el error en función del número de iteraciones para llegar a la convergencia en una semilogarítmica
axes[1].semilogy(errores, lw=2)
axes[1].set_xlabel('Iteración')
axes[1].set_ylabel('Error máximo [V]')
axes[1].set_title('Convergencia')
plt.tight_layout()
plt.savefig('V_r_convergencia.png', dpi=150)
print("Figura 1 guardada.")

# Verificamos que las CC se cumplen
print("V max:", np.max(V))
print("V min:", np.min(V))

# Mapa 2D
theta = np.linspace(0, 2*np.pi, 360) # ángulo polar
R2, TH = np.meshgrid(r, theta) # malla
X2 = R2 * np.cos(TH) # coordenada x
Y2 = R2 * np.sin(TH) # coordenada y
V2 = np.tile(V, (len(theta), 1)) # repetición angular (campo radial)

# niveles consistentes que van desde V_máx a V_mín para evitar que las interpolaciones tengan relevancia
levels = np.linspace(V.min(), V.max(), 40)

fig2, ax2 = plt.subplots(figsize=(7, 6.5))
cp = ax2.contourf(
    X2*100, Y2*100, V2,
    levels=levels,
    cmap='plasma'
)
# barra de color
cbar = plt.colorbar(cp, ax=ax2)

# valores "bonitos" de la barra de color: mínimo, máximo y puntos intermedios
cbar.set_ticks([
    V.min(),
    0.25 * V.max(),
    0.50 * V.max(),
    0.75 * V.max(),
    V.max()
])
cbar.set_label('V [V]')
ax2.set_aspect('equal')
ax2.set_xlabel('x [cm]')
ax2.set_ylabel('y [cm]')
ax2.set_title('Equipotenciales — sección transversal')
plt.tight_layout()
plt.savefig('equipotenciales_2D.png', dpi=150)
print("Fig 2 guardada.")