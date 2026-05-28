# Código que resuelve la práctica 10-11 y produce las simulaciones en formato
# .gif para aumentar la rapidez de ejecución

import numpy as np                                      # Operaciones numéricas y arrays
import matplotlib.pyplot as plt                         # Visualización de resultados
from matplotlib.animation import FuncAnimation, PillowWriter  # Animaciones y exportación a GIF
from scipy.linalg import solve_banded                   # Resolución eficiente de sistemas tridiagonales

# PARÁMETROS GLOBALES

L   = 1.0                          # Longitud de la caja [0, L]
dx  = 5e-4                         # Paso espacial de la discretización
x   = np.arange(0, L + dx, dx)    # Vector de posiciones: N puntos equiespaciados en [0, L]
Nx  = len(x)                       # Número total de puntos espaciales

dt  = 5e-6                         # Paso temporal de integración
Nt  = 5000                         # Número total de pasos temporales

x0    = 0.20                       # Posición inicial del centro del paquete gaussiano
sigma = 0.025                      # Anchura (dispersión) del paquete gaussiano
k0    = 120                        # Momento inicial del paquete (determina la velocidad media)
E     = k0**2 / 2                  # Energía cinética del paquete: E = k₀²/2 (unidades naturales, ℏ=m=1)

xc      = 0.60                     # Posición central de la barrera de potencial
V0_base = 1.2 * E                  # Altura de la barrera base: V0 > E → régimen de efecto túnel
w_base  = 0.015                    # Anchura de la barrera base

# FUNCIONES AUXILIARES

def construir_potencial(V0, w):
    """Barrera rectangular de altura V0 y anchura w centrada en xc."""
    V = np.zeros(Nx)                              # Potencial nulo en todo el dominio por defecto
    mask = (x >= xc - w/2) & (x <= xc + w/2)    # Máscara booleana: puntos dentro de la barrera
    V[mask] = V0                                  # Asigna V0 a los puntos interiores de la barrera
    return V


def inicializar_psi():
    """Paquete gaussiano normalizado con momento k0."""
    psi = np.exp(-(x - x0)**2 / (2 * sigma**2)) * np.exp(1j * k0 * x)
    # Gaussiana centrada en x0 con anchura sigma, modulada por onda plana e^{ik0·x}
    psi /= np.sqrt(np.trapezoid(np.abs(psi)**2, x))  # Normalización: ∫|ψ|²dx = 1
    psi[0]  = 0                                        # Condición de contorno: ψ(0) = 0 (pared izquierda)
    psi[-1] = 0                                        # Condición de contorno: ψ(L) = 0 (pared derecha)
    return psi


def construir_matrices_CN(V):
    """
    Construye las matrices A y B del esquema de Crank-Nicolson.

    El esquema resuelve: A·psi^n+1 = B·ψ^n
    A es tridiagonal implícita (avanzada en tiempo),
    B es tridiagonal explícita (conocida en el paso actual).

    Devuelve
    --------
    ab       : matriz A en formato bandeado (3 × Nx) para solve_banded
    lower_B  : subdiagonal de B  (longitud Nx-1)
    diag_B   : diagonal principal de B  (longitud Nx)
    upper_B  : superdiagonal de B  (longitud Nx-1)
    """
    alpha = 1j * dt / (4 * dx**2)     # Coeficiente complejo de la discretización temporal-espacial

    # Inicialización de diagonales como vectores de unos (los extremos quedan a 1 por las CC)
    diag_A  = np.ones(Nx,    dtype=complex)   # Diagonal principal de A
    diag_B  = np.ones(Nx,    dtype=complex)   # Diagonal principal de B
    upper_A = np.zeros(Nx-1, dtype=complex)   # Superdiagonal de A
    lower_A = np.zeros(Nx-1, dtype=complex)   # Subdiagonal de A
    upper_B = np.zeros(Nx-1, dtype=complex)   # Superdiagonal de B
    lower_B = np.zeros(Nx-1, dtype=complex)   # Subdiagonal de B

    for i in range(1, Nx-1):          # Recorre los puntos interiores (los extremos son CC fijas)

        # Elemento diagonal de A: término cinético (2alpha) + término potencial (idt·V_i/2)
        diag_A[i]    =  1 + 2*alpha + 1j*dt*V[i]/2

        # Elemento diagonal de B: igual que A pero con signos opuestos (esquema centrado en t)
        diag_B[i]    =  1 - 2*alpha - 1j*dt*V[i]/2

        # Elementos fuera de la diagonal (acoplamiento con vecinos): derivada segunda en diferencias finitas
        upper_A[i]   = -alpha          # A: acoplamiento con el punto i+1
        lower_A[i-1] = -alpha          # A: acoplamiento con el punto i-1
        upper_B[i]   =  alpha          # B: acoplamiento con el punto i+1 (signo opuesto a A)
        lower_B[i-1] =  alpha          # B: acoplamiento con el punto i-1 (signo opuesto a A)

    # Empaqueta A en formato bandeado requerido por solve_banded (3 filas: super, diag, sub)
    ab = np.zeros((3, Nx), dtype=complex)
    ab[0, 1:]  = upper_A               # Fila 0: superdiagonal (desplazada 1 posición a la derecha)
    ab[1, :]   = diag_A                # Fila 1: diagonal principal
    ab[2, :-1] = lower_A              # Fila 2: subdiagonal (desplazada 1 posición a la izquierda)

    return ab, lower_B, diag_B, upper_B


def apply_B(psi, lower_B, diag_B, upper_B):
    """
    Calcula el producto B·ψ sin construir la matriz B completa.
    Aprovecha la estructura tridiagonal para hacerlo en O(N).
    """
    result = np.zeros_like(psi, dtype=complex)     # Vector resultado inicializado a cero
    result[1:-1] = (
        lower_B[:-1] * psi[:-2]        # Contribución del vecino izquierdo  (subdiagonal)
        + diag_B[1:-1] * psi[1:-1]    # Contribución del punto propio       (diagonal)
        + upper_B[1:]  * psi[2:]       # Contribución del vecino derecho     (superdiagonal)
    )
    return result


def T_analitica(E, V0, w):
    """
    Transmisión exacta para una barrera rectangular.
    Usa la fórmula con sinh si E < V0 (túnel) o sin si E > V0 (sobre-barrera).
    """
    if E < V0:
        kappa = np.sqrt(2 * (V0 - E))             # Longitud de decaimiento evanescente dentro de la barrera
        return 1 / (1 + (V0**2 * np.sinh(kappa * w)**2) / (4 * E * (V0 - E)))
    else:
        q = np.sqrt(2 * (E - V0))                  # Momento real dentro de la barrera (E > V0)
        return 1 / (1 + (V0**2 * np.sin(q * w)**2)    / (4 * E * (E - V0)))


def evolucion_completa(V0, w, guardar_frames=False):
    """
    Ejecuta la evolución temporal completa con Crank-Nicolson.

    Parámetros
    ----------
    V0             : altura de la barrera de potencial
    w              : anchura de la barrera de potencial
    guardar_frames : si True, almacena frames de |psi|^2 y Re(psi) cada 20 pasos (para animación)

    Devuelve
    --------
    times      : instantes de tiempo en los que se muestreó T, R y P
    T_list     : probabilidad de transmisión en cada instante muestreado
    R_list     : probabilidad de reflexión en cada instante muestreado
    P_list     : probabilidad en la barrera en cada instante muestreado
    frames_prob: lista de arrays |psi|^2 (vacía si guardar_frames=False)
    frames_real: lista de arrays Re(psi) (vacía si guardar_frames=False)
    """
    V   = construir_potencial(V0, w)       # Genera el perfil de potencial con los parámetros dados
    psi = inicializar_psi()                # Estado cuántico inicial: paquete gaussiano normalizado
    ab, lB, dB, uB = construir_matrices_CN(V)  # Matrices del esquema CN para este potencial

    T_list, R_list, P_list = [], [], []    # Listas para almacenar T, R y P_barrera a lo largo del tiempo
    frames_prob, frames_real = [], []      # Listas para frames de animación (solo si se solicitan)
    times = []                             # Instantes de muestreo

    # Máscaras espaciales fijas para integrar en cada región
    left    = x < (xc - w/2)                               # Región izquierda: antes de la barrera
    barrier = (x >= xc - w/2) & (x <= xc + w/2)           # Región de la barrera
    right   = x > (xc + w/2)                               # Región derecha: después de la barrera

    for n in range(Nt):                    # Bucle principal de integración temporal

        b   = apply_B(psi, lB, dB, uB)    # Calcula el lado derecho: b = B·psi_n
        psi = solve_banded((1, 1), ab, b) # Resuelve A·psi_n+1 = b → avanza un paso temporal en O(N)
        psi[0] = psi[-1] = 0             # Reimpone condiciones de contorno (paredes impenetrables)
        psi /= np.sqrt(np.trapezoid(np.abs(psi)**2, x))  # Renormaliza para corregir deriva numérica acumulada

        if n % 20 == 0:                    # Muestrea cada 20 pasos para reducir el volumen de datos

            prob = np.abs(psi)**2                                        # Densidad de probabilidad

            T_list.append(np.trapezoid(prob[right],   x[right]))        # T = integral a la derecha de la barrera
            R_list.append(np.trapezoid(prob[left],    x[left]))         # R = integral a la izquierda de la barrera
            P_list.append(np.trapezoid(prob[barrier], x[barrier]))      # P = integral en la barrera
            times.append(n * dt)                                         # Guarda el instante de tiempo

            if guardar_frames:                     # Solo si se va a generar animación
                frames_prob.append(prob.copy())    # Guarda copia de |ψ|² para el frame
                frames_real.append(np.real(psi).copy())  # Guarda copia de Re(ψ) para el frame

    return times, T_list, R_list, P_list, frames_prob, frames_real


# Simulación base

print("Simulación base")

V_base = construir_potencial(V0_base, w_base)   # Potencial de la simulación base (para los plots)

# Ejecuta la evolución completa guardando frames para las animaciones
times, T_list, R_list, P_list, frames_prob, frames_real = evolucion_completa(
    V0_base, w_base, guardar_frames=True
)

T_arr = np.array(T_list)   # Convierte a array de NumPy para operar vectorialmente
R_arr = np.array(R_list)
P_arr = np.array(P_list)

# ---------- Resultados numéricos ----------
print(f"\nTransmisión final   T  = {T_arr[-1]:.6f}")                         # Último valor de T
print(f"Reflexión final     R  = {R_arr[-1]:.6f}")                           # Último valor de R
print(f"T + R               = {T_arr[-1]+R_arr[-1]:.6f}")                    # Debe ser ≈ 1 al final
print(f"T + R + P_barrera   = {T_arr[-1]+R_arr[-1]+P_arr[-1]:.6f}")         # Conservación total

# ---------- Resultado analítico ----------
T_teo = T_analitica(E, V0_base, w_base)   # Transmisión teórica exacta para estos parámetros
R_teo = 1 - T_teo                          # Reflexión teórica (R = 1 - T para barrera sin absorción)

print(f"\nEnergía       E  = {E:.3f}")
print(f"Potencial    V0  = {V0_base:.3f}")
print(f"Teoría        T  = {T_teo:.6f}")
print(f"Teoría        R  = {R_teo:.6f}")

# ---------- Comparación con máximo numérico ----------
idx_max = np.argmax(T_arr)         # Índice del instante en que T es máxima (paquete ya ha cruzado)
T_max   = T_arr[idx_max]           # Valor máximo de T durante la simulación

print(f"\nPrimer máximo numérico de T:")
print(f"  T_max            = {T_max:.6f}")
print(f"  R(T_max)         = {R_arr[idx_max]:.6f}")
print(f"  P_barrera(T_max) = {P_arr[idx_max]:.6f}")
print(f"  T+R+P            = {T_max+R_arr[idx_max]+P_arr[idx_max]:.6f}")
print(f"  t(T_max)         = {times[idx_max]:.4e}")
print(f"  T_teórica        = {T_teo:.6f}")
print(f"  Error relativo   = {abs(T_max-T_teo)/T_teo:.4%}")   # Discrepancia numérico vs analítico

# ============================================================
# Animaciones
# ============================================================

print("\nGenerando animaciones…")

# --- Animación de la densidad de probabilidad |ψ|² ---
fig1, ax1 = plt.subplots(figsize=(10, 5))
line1, = ax1.plot([], [], lw=2)                                   # Línea vacía que se actualizará
ax1.fill_between(x, 0, V_base / np.max(V_base) * 0.5,
                 color='red', alpha=0.3)                           # Perfil de la barrera escalado visualmente
txt1 = ax1.text(0.02, 0.90, '', transform=ax1.transAxes,
                fontsize=12)                                        # Cuadro de texto para T, R, P en tiempo real
ax1.set(xlim=(0, 1), ylim=(0, np.max(frames_prob) * 1.2),
        xlabel='x', ylabel=r'$|\psi(x,t)|^2$',
        title='Densidad de probabilidad')
ax1.grid(True)

def update1(f):
    """Función de actualización para cada frame de la animación |ψ|²."""
    line1.set_data(x, frames_prob[f])                              # Actualiza la curva con el frame f
    txt1.set_text(f"t = {times[f]:.5e}\nT = {T_list[f]:.4f}"
                  f"\nR = {R_list[f]:.4f}\nP_barrera = {P_list[f]:.4f}")
    return line1, txt1

ani1 = FuncAnimation(fig1, update1, frames=len(frames_prob),
                     interval=30, blit=False)                      # Genera la animación
ani1.save("densidad_probabilidad.gif", writer=PillowWriter(fps=25)) # Exporta a GIF a 25 fps
plt.close(fig1)                                                     # Cierra la figura para liberar memoria

# --- Animación de la parte real Re(psi) ---
fig2, ax2 = plt.subplots(figsize=(10, 5))
line2, = ax2.plot([], [], lw=2)
ax2.fill_between(x, -1, V_base / np.max(V_base),
                 color='red', alpha=0.3)                           # Barrera superpuesta al plot de Re(ψ)
txt2 = ax2.text(0.02, 0.90, '', transform=ax2.transAxes, fontsize=12)
ax2.set(xlim=(0, 1), ylim=(-1.2, 1.2),
        xlabel='x', ylabel=r'$\mathrm{Re}(\psi)$',
        title='Parte real de la función de onda')
ax2.grid(True)

def update2(f):
    """Función de actualización para cada frame de la animación Re(ψ)."""
    line2.set_data(x, frames_real[f])                              # Actualiza con la parte real del frame f
    txt2.set_text(f"t = {times[f]:.5e}")
    return line2, txt2

ani2 = FuncAnimation(fig2, update2, frames=len(frames_real),
                     interval=30, blit=False)
ani2.save("parte_real.gif", writer=PillowWriter(fps=25))
plt.close(fig2)

print("Animaciones guardadas: densidad_probabilidad.gif  |  parte_real.gif")

# ============================================================
# T(t), R(t), P_barrera(t) — conservación de probabilidad
# ============================================================

plt.figure(figsize=(9, 5))
plt.plot(times, T_arr,             label='T(t)',           lw=2)   # Transmisión acumulada en la región derecha
plt.plot(times, R_arr,             label='R(t)',           lw=2)   # Reflexión acumulada en la región izquierda
plt.plot(times, P_arr,             label='P_barrera(t)',   lw=2)   # Probabilidad transitoria dentro de la barrera
plt.plot(times, T_arr+R_arr+P_arr, '--', label='T+R+P_barrera', lw=2)  # Suma total: debe mantenerse ≈ 1
plt.xlabel('Tiempo')
plt.ylabel('Probabilidad')
plt.title('Conservación de probabilidad')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('conservacion_probabilidad.pdf', bbox_inches='tight')  # Guarda en PDF para el informe
plt.show()

# ============================================================
# Barrido en V0
# ============================================================

print("\n" + "=" * 50)
print("Barrido en V0")
print("=" * 50)

V0_values = V0_base * np.linspace(0.9, 5, 9)
T_map_V0  = []                                      # Lista para almacenar T(t) de cada simulación

for V0 in V0_values:
    print(f"   V0 = {V0:.1f}…", end='', flush=True)
    _, T_v, _, _, _, _ = evolucion_completa(V0, w_base)   # Simula con este V0 y w fijo
    T_map_V0.append(T_v)                                    # Acumula la serie T(t) resultante
    print(" done")

T_map_V0 = np.array(T_map_V0)   # Convierte a matriz 2D: filas = valores de V0, columnas = instantes

# --- Mapa de calor T(V0, t) ---
plt.figure(figsize=(8, 5))
plt.imshow(T_map_V0, aspect='auto', origin='lower',
           extent=[0, Nt*dt, V0_values[0], V0_values[-1]])  # Ejes: tiempo en x, V0 en y
plt.colorbar(label='T(t)')
plt.xlabel('Tiempo')
plt.ylabel('$V_0$')
plt.title('Evolución $T(V_0, t)$')
plt.tight_layout()
plt.savefig('mapa_T_V0.pdf', bbox_inches='tight')
plt.show()

# --- Promedio temporal ± desviación estándar ---
T_mean_V0 = T_map_V0.mean(axis=1)   # Media temporal de T para cada valor de V0
T_std_V0  = T_map_V0.std(axis=1)    # Desviación estándar temporal (barra de error)

plt.figure(figsize=(7, 5))
plt.errorbar(V0_values, T_mean_V0, yerr=T_std_V0,
             marker='o', capsize=4, linestyle='-')    # Barras de error = variabilidad temporal de T
plt.xlabel('$V_0$')
plt.ylabel(r'$\langle T \rangle$ (promedio temporal)')
plt.title('Transmisión media vs altura de barrera')
plt.grid(True)
plt.tight_layout()
plt.savefig('T_medio_vs_V0.pdf', bbox_inches='tight')
plt.show()

# --- Sin barras de error ---
plt.figure(figsize=(7, 5))
plt.plot(V0_values, T_mean_V0, marker='o')            # Misma gráfica pero más limpia visualmente
plt.xlabel('$V_0$')
plt.ylabel(r'$\langle T \rangle$ (promedio temporal)')
plt.title('Transmisión media vs altura de barrera (sin σ)')
plt.grid(True)
plt.tight_layout()
plt.savefig('T_medio_vs_V0_noerror.pdf', bbox_inches='tight')
plt.show()

# ============================================================
# Barrido en w
# ============================================================

print("\n" + "=" * 50)
print("Barrido en w")
print("=" * 50)

w_values = w_base * np.linspace(0.2, 5, 9)   # 9 valores de w entre 0.5·w_base y 1.5·w_base
T_map_w  = []                                    # Lista para almacenar T(t) de cada simulación

for w in w_values:
    print(f"   w = {w:.5f}…", end='', flush=True)
    _, T_w, _, _, _, _ = evolucion_completa(V0_base, w)   # Simula con V0 fijo y este valor de w
    T_map_w.append(T_w)                                     # Acumula la serie T(t) resultante
    print(" done")

T_map_w = np.array(T_map_w)   # Convierte a matriz 2D: filas = valores de w, columnas = instantes

# --- Mapa de calor T(w, t) ---
plt.figure(figsize=(8, 5))
plt.imshow(T_map_w, aspect='auto', origin='lower',
           extent=[0, Nt*dt, w_values[0], w_values[-1]])   # Ejes: tiempo en x, w en y
plt.colorbar(label='T(t)')
plt.xlabel('Tiempo')
plt.ylabel('$w$')
plt.title('Evolución $T(w, t)$')
plt.tight_layout()
plt.savefig('mapa_T_w.pdf', bbox_inches='tight')
plt.show()

# --- Promedio temporal ± desviación estándar ---
T_mean_w = T_map_w.mean(axis=1)   # Media temporal de T para cada valor de w
T_std_w  = T_map_w.std(axis=1)    # Desviación estándar temporal

plt.figure(figsize=(7, 5))
plt.errorbar(w_values, T_mean_w, yerr=T_std_w,
             marker='o', capsize=4, linestyle='-')    # La fuerte caída exponencial con w es la firma del túnel
plt.xlabel('$w$')
plt.ylabel(r'$\langle T \rangle$ (promedio temporal)')
plt.title('Transmisión media vs anchura de barrera')
plt.grid(True)
plt.tight_layout()
plt.savefig('T_medio_vs_w.pdf', bbox_inches='tight')
plt.show()

# --- Sin barras de error ---
plt.figure(figsize=(7, 5))
plt.plot(w_values, T_mean_w, marker='o')
plt.xlabel('$w$')
plt.ylabel(r'$\langle T \rangle$ (promedio temporal)')
plt.title('Transmisión media vs anchura de barrera (sin σ)')
plt.grid(True)
plt.tight_layout()
plt.savefig('T_medio_vs_w_noerror.pdf', bbox_inches='tight')
plt.show()

print("\nTodos los resultados generados.")