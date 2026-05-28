import numpy as np
import matplotlib.pyplot as plt
from time import time

# no se puede usar numba con time y entonces el código es muy lento

# Parámetros
M      = 100     # grid squares per side
V      = 1.0     # top-wall voltage [V]
target = 1e-6    # convergence criterion [V]


def solve_laplace_sor(omega, target=1e-6, M=100, V=1.0):
    """
    Solve Laplace's equation with SOR parameter omega.
    Returns (phi, iterations, elapsed_time).
    """
    phi = np.zeros((M + 1, M + 1), float)
    phi[M, :] = V     # top wall  (row M = top when origin is bottom)

    delta      = 1.0
    iterations = 0
    t0         = time()

    coeff_sum  = (1.0 + omega) / 4.0
    coeff_self = omega

    while delta > target:
        delta = 0.0
        # Gauss-Seidel usa los valores nuevos instantáneamente
        for i in range(1, M):
            for j in range(1, M):
                new = coeff_sum * (
                    phi[i + 1, j] + phi[i - 1, j] +
                    phi[i, j + 1] + phi[i, j - 1]
                ) - coeff_self * phi[i, j]

                diff = abs(new - phi[i, j])
                if diff > delta:
                    delta = diff
                phi[i, j] = new

        iterations += 1

    return phi, iterations, time() - t0


# Puro Gauss-Seidel
print("Resolviendo w=0 …")
phi_gs, it_gs, t_gs = solve_laplace_sor(omega=0.0)
print(f"  Iterations: {it_gs}   Time: {t_gs:.1f} s") # 200.7 s

# Sobre-relajación cerca del valor óptimo de omega
print("\nResolviendo w=0.9 …")
phi_sor, it_sor, t_sor = solve_laplace_sor(omega=0.9)
print(f"  Iterations: {it_sor}   Time: {t_sor:.1f} s")
print(f"  Speed-up over ω=0: ×{it_gs / it_sor:.1f}") # 15.9 s

# Buscar w óptimo
print("\nEscaneando w …")
omega_values = np.arange(0.5, 1.0, 0.1)
iter_counts  = []

for w in omega_values:
    _, iters, _ = solve_laplace_sor(omega=w)
    iter_counts.append(iters)
    print(f"  ω = {w:.2f}  →  {iters} iterations")

omega_opt  = omega_values[np.argmin(iter_counts)]
iters_opt  = min(iter_counts)
print(f"\n  w óptimo = {omega_opt:.2f}  ({iters_opt} iterations)")

# Re-solve with optimal ω for the final plot
phi_opt, _, _ = solve_laplace_sor(omega=omega_opt)

# ── Plotting ─────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(15, 5))
fig.suptitle(
    "Laplace's equation – Exercise 9.2  (Gauss–Seidel + over-relaxation)",
    fontsize=12, fontweight='bold'
)

extent = [0, M, M, 0]   # indices, origin top-left (matches imshow default)

# ---- Panel 1: φ with ω = 0 (Gauss-Seidel) ----------------------------------
ax1 = fig.add_subplot(1, 3, 1)
im1 = ax1.imshow(phi_gs, cmap='jet', extent=extent)
plt.colorbar(im1, ax=ax1, label='φ (V)')
ax1.set_title(f'Pure Gauss–Seidel  (ω = 0)\n{it_gs} iterations', fontsize=10)
ax1.set_xlabel('x  (grid units)')
ax1.set_ylabel('y  (grid units)')

# ---- Panel 2: φ with optimal ω (SOR) ---------------------------------------
ax2 = fig.add_subplot(1, 3, 2)
im2 = ax2.imshow(phi_opt, cmap='jet', extent=extent)
plt.colorbar(im2, ax=ax2, label='φ (V)')
ax2.set_title(
    f'SOR  (ω = {omega_opt:.2f}, optimal)\n{iters_opt} iterations', fontsize=10
)
ax2.set_xlabel('x  (grid units)')

# ---- Panel 3: Iteration count vs ω -----------------------------------------
ax3 = fig.add_subplot(1, 3, 3)
ax3.plot(omega_values, iter_counts, 'o-', color='steelblue',
         markersize=6, linewidth=2, label='Iterations')
ax3.axvline(omega_opt, color='red', linestyle='--',
            label=f'Optimal ω = {omega_opt:.2f}')
ax3.fill_between(omega_values, iter_counts, alpha=0.15, color='steelblue')
ax3.set_xlabel('Relaxation parameter ω', fontsize=11)
ax3.set_ylabel('Iterations to convergence', fontsize=11)
ax3.set_title('Convergence speed vs ω\n(lower = faster)', fontsize=10)
ax3.legend(fontsize=10)
ax3.grid(True, alpha=0.4)

plt.tight_layout()
plt.savefig('/mnt/user-data/outputs/ejercicio_9_2.png', dpi=150, bbox_inches='tight')
plt.show()
print("\nFigure saved → ejercicio_9_2.png")
