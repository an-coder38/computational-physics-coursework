# condensador con cilindro interno a 2.5V y externo a tierra, las tapas en 0 y L están a 0V -- problema depende de z, independiente de theta
import numpy as np
import matplotlib.pyplot as plt
from numba import njit

# Parámetros
Ri, Re = 0.01, 0.04
L = 0.20
Vi = 2.5

# Malla 2D (r,z)
Nr, Nz = 120, 120

r = np.linspace(Ri, Re, Nr)
z = np.linspace(0, L, Nz)

dr = r[1] - r[0]
dz = z[1] - z[0]

# Potencial inicializado
V = np.zeros((Nr, Nz))

# condición inicial: solo conductor interno
for i in range(Nr):
    V[i, :] = Vi * (Re - r[i]) / (Re - Ri)

epsilon = 1e-6
max_iter = 50000


# Gauss-Seidel 2D en cilíndricas
@njit
def solve(V, r, dr, dz, Ri, Re, Vi, epsilon, max_iter):

    Nr, Nz = V.shape
    errores = np.zeros(max_iter)

    for it in range(max_iter):

        delta = 0.0

        for i in range(1, Nr - 1):
            for j in range(1, Nz - 1):

                old = V[i, j]

                ri = r[i]

                # término radial (con geometría cilíndrica)
                Vr_plus  = V[i + 1, j]
                Vr_minus = V[i - 1, j]

                fac = dr / (2.0 * ri)
                radial = (1 - fac) * Vr_minus + (1 + fac) * Vr_plus
                radial *= 0.5

                # término axial
                Vz = (V[i, j + 1] + V[i, j - 1]) / 2.0

                # combinación Laplace 2D
                V[i, j] = 0.5 * (radial + Vz)

                diff = abs(V[i, j] - old)
                if diff > delta:
                    delta = diff

        # CC
        V[0, :] = Vi # interno
        V[-1, :] = 0.0 #externo
        V[:, 0] = 0.0 # tapa
        V[:, -1] = 0.0 # tapa

        errores[it] = delta

        if delta < epsilon:
            return V, errores[:it+1], it+1

    return V, errores, max_iter


# Usar la función
V, errores, it = solve(V, r, dr, dz, Ri, Re, Vi, epsilon, max_iter)

print("Convergencia en:", it)
print("V max:", V.max(), "V min:", V.min())


# Mapa V(z,r)
R, Z = np.meshgrid(r, z, indexing='ij')
fig, ax = plt.subplots(figsize=(7, 6))
levels = np.linspace(0, Vi, 40)
cp = ax.contourf(
    R*100, Z*100, V,
    levels=levels,
    cmap='plasma'
)
cbar = plt.colorbar(cp, ax=ax)
cbar.set_ticks([0,0.25*Vi, 0.5*Vi,0.75*Vi, Vi])
ax.set_xlabel('r [cm]')
ax.set_ylabel('z [cm]')
cbar.set_label('V [V]')
ax.set_title('Potencial 2D con tapas a tierra')
plt.tight_layout()
plt.savefig('equipotenciales_2D_tierra.png', dpi=150)
plt.show()