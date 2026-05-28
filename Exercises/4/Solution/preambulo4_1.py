# Llamada a las librerías y funciones necesarias
import math as m
import numpy as np
from pylab import plot, show, xlabel, ylabel, title, legend


def pendulo_no_lineal(q, FI, OmegaI):

    theta = np.zeros(N)
    omega = np.zeros(N)
    t = np.zeros(N)

    theta[0] = theta0
    omega[0] = omega0
    t[0] = t0

    for i in range(N-1):

        a = - (g/l) * np.sin(theta[i]) \
            - q * omega[i] \
            + FI * np.sin(OmegaI * t[i])

        omega[i+1] = omega[i] + a * dt
        theta[i+1] = theta[i] + omega[i+1] * dt
        t[i+1] = t[i] + dt

    return t, theta


# Inicialización de variables
g = 9.8
l = 1
theta0 = m.radians(5) # Permite la entrada del ángulo inicial en grados
omega0 = 0.
t0 = 0.

timeTotal = 30

dt = 0.04           # Intervalo de tiempos
NN = timeTotal / dt # Número de elementos de las colecciones
N = int(NN)         # Transformación de NN en un entero.

theta = np.zeros(N, float) # Declaración de las colecciones a utilizar: ángulo
omega = np.zeros(N, float) # Velocidad angular
t = np.zeros(N,float)      # Tiempo

# Valores iniciales de las colecciones
theta[0] = theta0
omega[0] = omega0
t[0] = t0

# Bucle que calcula los valores de las coordenadas "i+1" en función de las de "i".
for i in range(N-1):
    omega[i+1] = omega[i] - (g/l)*theta[i]*dt
    theta[i+1] = theta[i] + omega[i+1]*dt           #Método de Euler-Cromer
    t[i+1] = t[i] + dt

plot(t,theta)
xlabel("t (s)"), ylabel(r'$\theta$ (rad)')
title("Péndulo Ideal - Método de Euler-Cromer")
show()

q = 0.5 # amortiguado 4, sobreamortiguado 10
theta[0] = theta0
omega[0] = omega0
t[0] = t0
for i in range(N-1):
    omega[i+1] = omega[i] - (g/l)*theta[i]*dt - q * omega[i]*dt
    theta[i+1] = theta[i] + omega[i+1]*dt           #Método de Euler-Cromer
    t[i+1] = t[i] + dt

# Instrucciones para la realización de las gráficas.
# an = sin(OmegaI*t+phi)
t_no, theta_no = pendulo_no_lineal(q, 0, 0)
plot(t,theta, label = 'lineal')
plot(t_no, theta_no, label='no lineal')
legend()
xlabel("t (s)"), ylabel(r'$\theta$ (rad)')
title("Péndulo subamortiguado - Método de Euler-Cromer")
show()

q = 0.5 # amortiguado 4, sobreamortiguado 10
OmegaI=0.5
phi= np.radians(30)
FI = 0.01
theta[0] = theta0
omega[0] = omega0
t[0] = t0
for i in range(N-1):
    omega[i+1] = omega[i] - (g/l)*theta[i]*dt - q * omega[i]*dt + FI * np.sin(OmegaI*t[i])*dt
    theta[i+1] = theta[i] + omega[i+1]*dt           #Método de Euler-Cromer
    t[i+1] = t[i] + dt

# Instrucciones para la realización de las gráficas.
# an = sin(OmegaI*t+phi)
t_no, theta_no = pendulo_no_lineal(q, OmegaI, FI)
plot(t,theta, label = 'lineal')
plot(t_no, theta_no, label='no lineal')
legend()
xlabel("t (s)"), ylabel(r'$\theta$ (rad)')
title("Con fuerza impulsora")
show()

q = 0.5 # amortiguado 4, sobreamortiguado 10
OmegaI = np.sqrt(g/l)
phi= np.radians(30)
FI = 0.01
theta[0] = theta0
omega[0] = omega0
t[0] = t0
for i in range(N-1):
    omega[i+1] = omega[i] - (g/l)*theta[i]*dt - q * omega[i]*dt + FI * np.sin(OmegaI*t[i])
    theta[i+1] = theta[i] + omega[i+1]*dt           #Método de Euler-Cromer
    t[i+1] = t[i] + dt

# Instrucciones para la realización de las gráficas.
# an = sin(OmegaI*t+phi)
plot(t,theta, label = 'lineal')
legend()
xlabel("t (s)"), ylabel(r'$\theta$ (rad)')
title("Con fuerza impulsora cerca frecuencia natural")
show()
