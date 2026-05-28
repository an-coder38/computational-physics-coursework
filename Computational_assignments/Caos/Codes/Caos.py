import numpy as np
import matplotlib.pyplot as plt


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

        a = -(g/l)*np.sin(theta[i]) - q*omega[i] + FI*np.sin(OmegaI*t[i])

        omega[i+1] = omega[i] + a*dt
        theta[i+1] = theta[i] + omega[i+1]*dt
        t[i+1] = t[i] + dt

    return t, theta

def delta_pendulo_no_lineal (theta0, dtheta, omega0, timeTotal, t0, dt, q, FI, OmegaI):
    # calcula el valor absoluto de un ángulo menos uno cercano (cercanía dtheta) en función del tiempo
    
    t, theta = pendulo_no_lineal(theta0, omega0, timeTotal, t0, dt, q, FI, OmegaI)
    
    theta_mod=theta0+dtheta
    t, theta1 = pendulo_no_lineal(theta_mod, omega0, timeTotal, t0, dt, q, FI, OmegaI)
    delta = abs(theta-theta1)
    
    return t, delta

def delta_pendulo_no_lineal_q (theta0, dtheta, omega0, timeTotal, t0, dt, q, dq, FI, OmegaI):
    # calcula el valor absoluto de un ángulo menos uno cercano (cercanía dtheta) en función del tiempo
    
    t, theta = pendulo_no_lineal(theta0, omega0, timeTotal, t0, dt, q, FI, OmegaI)
    t, thetaq = pendulo_no_lineal(theta0+dtheta, omega0, timeTotal, t0, dt, q+dq, FI, OmegaI)
    delta = abs(theta-thetaq)
    
    return t, delta

def ajustes_lyapunov(t, delta):

    mask = delta > 1e-10
    t = t[mask]
    delta = delta[mask]

    mask_time = t > 15 #máscara para el tiempo de estado transitorio
    t = t[mask_time]
    delta = delta[mask_time]

    if len(t) < 2:   # evitar error en polyfit si no hay datos que cumplen ambas máscaras
        return t, np.array([]), np.nan, np.nan, np.array([])

    log_delta = np.log(delta)

    coef = np.polyfit(t, log_delta, 1)
    pendiente = coef[0]
    ordenada = coef[1]

    log_delta_fit = pendiente*t + ordenada
    return t, log_delta, pendiente, ordenada, log_delta_fit

# Condiciones iniciales
g = 9.8
l = 9.8
theta0 = 0.2
theta1 = 0.4
theta2 = 0.6
dtheta = 0.01
omega0 = 0.
t0 = 0.
timeTotal = 100
dt = 0.01           # Intervalo de tiempos


q_nula = 0
q_sub = 0.5 # amortiguado 4, sobreamortiguado 10
q_amo = 4.
q_sob = 10.
OmegaI= 2/3
FI = 1.2

pendientes=[]

# a) Exponente de Lyapunov con CI comparando theta y theta'=theta+dtheta
print('q fijo = ', q_sub)
for theta in [theta0, theta1, theta2]:
    t, delta = delta_pendulo_no_lineal(theta, dtheta, omega0, timeTotal, t0, dt, q_sub, FI, OmegaI)
    
    t, log_delta, pendiente, ordenada, log_delta_fit = ajustes_lyapunov (t,delta)
    print(f"Exponente de Lyapunov estimado para theta {theta:.2f} =", pendiente)

    #gráficas
    plt.figure(figsize=(8,5))
    plt.plot(t, log_delta, label="log(Δ(t))")
    plt.plot(t, log_delta_fit, '--', label=f"Ajuste lineal\npendiente = {pendiente:.4f}")
    plt.xlabel("t")
    plt.ylabel("log(Δ)")
    plt.legend()
    plt.grid()
    plt.savefig(f"apartado_a_theta_{theta:.2f}.png", dpi=300)
    plt.show()
    plt.close()
    
    
    #Ver el comportamiento de los péndulos al varial ligeramente el ángulo inicial
    ts, thetas = pendulo_no_lineal(theta, omega0, timeTotal, t0, dt, q_sub, FI, OmegaI) #theta0
    td, thetad = pendulo_no_lineal(theta+dtheta, omega0, timeTotal, t0, dt, q_sub, FI, OmegaI) #theta0+dtheta
    
    #gráficas
    plt.figure(figsize=(8,5))
    plt.plot(ts, thetas, label=f"theta = {theta:.2f}")
    plt.plot(td, thetad, label=f"theta = {theta+dtheta:.2f}")
    plt.xlabel("t (s)")
    plt.ylabel(r'$\theta$ (rad)')
    plt.legend()
    plt.grid()
    plt.savefig(f"comp_ci_{theta:.2f}.png", dpi=300)
    plt.show()
    plt.close()

    pendientes.append(pendiente)    

promedio = np.mean(pendientes)
desv = np.std(pendientes)

print(f"Exponente de Lyapunov estimado = {promedio:.5f} ± {desv:.5f}")




# b) Sensibilidad al amortiguamiento
# posiciones iniciales casi idénticas y pequeña variación en q
print('q variable')
dq = 0.01

for theta in [theta0, theta1, theta2]:
    t, delta = delta_pendulo_no_lineal_q(theta, dtheta, omega0, timeTotal, t0, dt, q_sub, dq, FI, OmegaI)
    
    t, log_delta, pendiente, ordenada, log_delta_fit = ajustes_lyapunov(t, delta)
    
    print("\nExponente de Lyapunov con diferencia en q =", pendiente)
    
    plt.figure(figsize=(8,5))
    
    plt.plot(t, log_delta, label="log(Δ(t))")
    plt.plot(t, log_delta_fit,'--', label="ajuste lineal")
    
    plt.xlabel("t")
    plt.ylabel("log(Δ)")
    plt.legend()
    plt.grid()
    plt.savefig(f"apartado_b_theta_{theta:.2f}.png", dpi=300)
    plt.show()


    #Ver el comportamiento de los péndulos al varial ligeramente el ángulo inicial
    ts, thetas = pendulo_no_lineal(theta, omega0, timeTotal, t0, dt, q_sub, FI, OmegaI) #theta0
    td, thetad = pendulo_no_lineal(theta+dtheta, omega0, timeTotal, t0, dt, q_sub+dq, FI, OmegaI) #theta0+dtheta
    
    #gráficas
    plt.figure(figsize=(8,5))
    plt.plot(ts, thetas, label=f"theta = {theta:.2f}")
    plt.plot(td, thetad, label=f"theta = {theta+dtheta:.2f}")
    plt.xlabel("t (s)")
    plt.ylabel(r'$\theta$ (rad)')
    plt.legend()
    plt.grid()
    plt.savefig(f"comp_ci_b_{theta:.2f}.png", dpi=300)
    plt.show()
    plt.close()

    pendientes.append(pendiente)    

promedio = np.mean(pendientes)
desv = np.std(pendientes)

print(f"Exponente de Lyapunov estimado = {promedio:.5f} ± {desv:.5f}")







# En este caso, se ve que el exponente de lyapunov se mantiene positivo, es decir, el sitema se mantiene en régimen caótico



# Extra: Estudio general
# Graficar evolución del exponente de lyapunov al aumentar q, para ver cuánto afecta el cambio en q
# y si su aumento o disminución puede cambiar el régimen caótico del sistema
# además de esta forma se aisla el comportamiento en función de q en vez de cambiar al mismo tiempo las
# condiciones iniciales y q
pendientes_list = []
q_super_sob = 100 # sistema extremadamente amortiguado
q_list=np.linspace(q_nula, q_super_sob,200) #lista con baja precisión pero un gran rango

for q in q_list:
    
    t, delta = delta_pendulo_no_lineal (theta0, dtheta, omega0, timeTotal, t0, dt, q, FI, OmegaI)
    
    t, log_delta, pendiente, ordenada, log_delta_fit = ajustes_lyapunov(t, delta)
    
    pendientes_list.append(pendiente)

plt.figure(figsize=(8,5))

plt.plot(q_list, pendientes_list)

plt.xlabel("q")
plt.ylabel("lambda")
plt.title("Evolución del exponente de lyapunov en función de q")
plt.grid()
plt.savefig("lyapunov(q).png", dpi=300)
plt.show()

# Zoom del exponente de lyapunov al aumentar q para ver qué está pasando sobre todo antes
# de que el sistema esté completamente dominado por q
pendientes_list = []
q_zoom = 3
q_list=np.linspace(q_nula, q_zoom,100)

for q in q_list:
    
    t, delta = delta_pendulo_no_lineal (theta0, dtheta, omega0, timeTotal, t0, dt, q, FI, OmegaI)
    
    t, log_delta, pendiente, ordenada, log_delta_fit = ajustes_lyapunov(t, delta)
    
    pendientes_list.append(pendiente)

plt.figure(figsize=(8,5))

plt.plot(q_list, pendientes_list)

plt.xlabel("q")
plt.ylabel("lambda")
plt.title("Evolución del exponente de lyapunov en función de q")
plt.grid()
plt.savefig("zoom_general_lyapunov(q).png", dpi=300)
plt.show()

# Zoom del exponente de lyapunov al aumentar q en la cercanía de 0.5
pendientes_list = []
q_sub_pre = 0.3
q_sub_post = 0.7
q_list=np.linspace(q_sub_pre, q_sub_post,100)

for q in q_list:
    
    t, delta = delta_pendulo_no_lineal (theta0, dtheta, omega0, timeTotal, t0, dt, q, FI, OmegaI)
    
    t, log_delta, pendiente, ordenada, log_delta_fit = ajustes_lyapunov(t, delta)
    
    pendientes_list.append(pendiente)

plt.figure(figsize=(8,5))

plt.plot(q_list, pendientes_list)

plt.xlabel("q")
plt.ylabel("lambda")
plt.title("Evolución del exponente de lyapunov en función de q")
plt.grid()
plt.savefig("zoom_qsub_lyapunov(q).png", dpi=300)
plt.show()
