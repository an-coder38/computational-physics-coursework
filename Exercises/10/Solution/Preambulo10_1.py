import numpy as np
import matplotlib.pyplot as plt

# Parámetros físicos
tau = 365.0 # Periodo orbital en días
A = 10.0 # Temperatura media superficial (°C)
B = 12.0 # Amplitud de variación térmica superficial (°C)
L = 20.0 # Profundidad máxima (m)
T_fondo = 11.0 # Temperatura a 20 m de profundidad (°C)
D = 0.1 # Difusividad térmica de la corteza (m^2 / día)

# Mallado
N = 100 # Número de divisiones espaciales
dx = L / N # Tamaño del paso espacial (0.2 m)
h = 0.1 # Tamaño del paso temporal (0.1 días)

# Constante de la ecuación discretizada
c = D * h / (dx**2) # c=0.25<0.5 que es necesario por Courante para que sea estable

# CI
T = np.full(N + 1, 10.0) # inicia en T salvo en el fondo
T[N] = T_fondo 
x = np.linspace(0, L, N + 1) # Array de profundidades

pasos_totales = int((10 * 365) / h) # Pasos para simular 10 años
pasos_9_anos = int((9 * 365) / h) # Momento donde inicia el 10mo año
pasos_trimestre = int((365 / 4) / h) # Pasos equivalentes a 3 meses

plt.figure(figsize=(10, 8))

for n in range(pasos_totales + 1):
    # simula para todos los años pero solo toma acción de guardar los datos en el décimo año
    t = n * h
    
    # Actualizar la condición de frontera en la superficie (x = 0)
    T[0] = A + B * np.sin(2 * np.pi * t / tau)
    
    # Extraer y graficar los perfiles durante el décimo año
    if n >= pasos_9_anos and (n - pasos_9_anos) % pasos_trimestre == 0: # solo que sea en los trimestres
        mes = int((n - pasos_9_anos) / pasos_trimestre) * 3
        if mes <= 9:  # Graficar a los 0, 3, 6 y 9 meses del último año
            plt.plot(x, T, label=f'T = {mes} meses (Año 10)', linewidth=2)
            
    # Calcular las nuevas temperaturas en el interior usando el esquema FTCS
    T_new = np.copy(T)
    T_new[1:N] = T[1:N] + c * (T[2:N+1] - 2*T[1:N] + T[0:N-1])
    
    # Actualizar el arreglo para el siguiente paso de tiempo
    T = np.copy(T_new)

plt.ylabel('Temperatura (°C)', fontsize=12)
plt.xlabel('Profundidad (m)', fontsize=12)
plt.title('Propagación del Calor en la Corteza Terrestre (Ciclo Anual)', fontsize=14)
plt.legend(loc='lower left', fontsize=10)
plt.show()