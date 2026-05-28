# Llamada a las librerías y funciones necesarias
import math as m
import numpy as np
import matplotlib.pyplot as plt

def pendulo_lineal(theta0, omega0, timeTotal, t0, dt, q, FI, OmegaI):
    #theta0 ya en radianes
    
    NN = (timeTotal-t0) / dt # Número de elementos de las colecciones
    N = int(NN)         # Transformación de NN en un entero.
    theta = np.zeros(N)
    omega = np.zeros(N)
    t = np.zeros(N)
    theta[0] = theta0
    omega[0] = omega0
    t[0] = t0


    for i in range(N-1):

        a = - (g/l) * theta[i] - q * omega[i] + FI * np.sin(OmegaI * t[i])

        omega[i+1] = omega[i] + a * dt
        theta[i+1] = theta[i] + omega[i+1] * dt
        t[i+1] = t[i] + dt

    return t, theta

def pendulo_no_lineal(theta0, omega0, timeTotal, t0, dt, q, FI, OmegaI):
    #theta0 ya en radianes
    
    NN = (timeTotal-t0) / dt # Número de elementos de las colecciones
    N = int(NN)         # Transformación de NN en un entero.
    theta = np.zeros(N)
    omega = np.zeros(N)
    t = np.zeros(N)
    theta[0] = theta0
    omega[0] = omega0
    t[0] = t0

    for i in range(N-1):

        a = - (g/l) * np.sin(theta[i]) - q * omega[i] + FI * np.sin(OmegaI * t[i])

        omega[i+1] = omega[i] + a * dt
        theta[i+1] = theta[i] + omega[i+1] * dt
        t[i+1] = t[i] + dt

    return t, theta

def delta_pendulo_no_lineal (theta0, dtheta, omega0, timeTotal, t0, dt, q, FI, OmegaI):
    # calcula el valor absoluto de un ángulo menos uno cercano (cercanía dtheta) en función del tiempo
    
    t, theta = pendulo_no_lineal(theta0, omega0, timeTotal, t0, dt, q, FI, OmegaI)
    theta_mod=theta0+dtheta
    t, theta1 = pendulo_no_lineal(theta_mod, omega0, timeTotal, t0, dt, q, FI, OmegaI)
    delta = abs(theta-theta1)
    
    return t, delta


# Condiciones iniciales
g = 9.8
l = 1
theta0 = 0.2 # Permite la entrada del ángulo inicial en grados
theta1 = 0.3
theta2 = 0.4
dtheta = 0.01
omega0 = 0.
t0 = 0.
timeTotal = 30
dt = 0.04           # Intervalo de tiempos


q_nula = 0
q_sub = 0.5 # amortiguado 4, sobreamortiguado 10
q_amo = 4
q_sob = 10
OmegaI= 2/3
FI = 4

pendientes=[]
for theta in [theta0, theta1, theta2]:
    t, delta = delta_pendulo_no_lineal(theta, dtheta, omega0, timeTotal, t0, dt, q_sub, FI, OmegaI)
    
    # Evitar log(0) con una máscara que elimine delta aprox igual a 0
    mask = delta > 1e-12
    t = t[mask]
    delta = delta[mask]
    
    log_delta = np.log(delta)
    
    coef = np.polyfit(t, log_delta, 1)
    
    pendiente = coef[0]
    ordenada = coef[1]
    
    print(f"Exponente de Lyapunov estimado para theta {theta} =", pendiente)
    
    # Recta ajustada
    log_delta_fit = pendiente * t + ordenada
    
    #gráficas
    plt.figure(figsize=(8,5))
    plt.plot(t, log_delta, label="log(Δ(t))")
    plt.plot(t, log_delta_fit, '--', label=f"Ajuste lineal\npendiente = {pendiente:.4f}")
    plt.xlabel("t")
    plt.ylabel("log(Δ)")
    plt.legend()
    plt.title(f'theta = {theta}')
    plt.grid()
    plt.show()
    plt.close()
    
    pendientes.append(pendiente)
    
promedio = np.mean(pendientes)
print(f"Exponente de Lyapunov promedio estimado =", promedio)



# Pendulo subamortiguado, amortiguado y sobreamortiguado
# for q in [q_nula, q_sub, q_amo, q_sob]:
#     for imp in ['sin', 'con']:
#         if imp=='sin':
#             t, theta = pendulo_lineal(theta0, omega0, timeTotal, t0, dt, q, 0, 0)
#             t_no, theta_no = pendulo_no_lineal(theta0, omega0, timeTotal, t0, dt, q, 0, 0)
#         if imp=='con':
#             t, theta = pendulo_lineal(theta0, omega0, timeTotal, t0, dt, q, FI, OmegaI)
#             t_no, theta_no = pendulo_no_lineal(theta0, omega0, timeTotal, t0, dt, q, FI, OmegaI)
#         plt.plot(t,theta, label = 'lineal')
#         plt.plot(t_no, theta_no, label='no lineal')
#         plt.legend()
#         plt.xlabel("t (s)"), plt.ylabel(r'$\theta$ (rad)')
#         if q==q_sub:
#            plt.title('Péndulo subamortiguado '+imp+' fuerza impulsora')
#         if q==q_amo:
#             plt.title('Péndulo amortiguado '+imp+' fuerza impulsora')    
#         if q==q_sob:
#             plt.title('Péndulo sobreamortiguado '+imp+' fuerza impulsora')    
#         plt.show()
#         if q==q_nula:
#             plt.title('Péndulo sin rozamiento '+imp+' fuerza impulsora')
    