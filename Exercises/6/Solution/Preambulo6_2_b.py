import numpy as np
import matplotlib.pyplot as plt
from random import random, seed

def fit_recta (x_min, x_max, m):
    x_fit = np.linspace(x_min, x_max, 20)
    y_fit = m*x_fit
    return x_fit, y_fit

m = 1000      # numero de caminantes
N = 500       # numero de pasos

np.random.seed(10)

pos = np.zeros(m)
msd = []

for n in range(N):
    
    pasos = np.random.choice([-1,1], size=m)
    pos += pasos
    
    msd.append(np.mean(pos**2))

x_fit, y_fit = fit_recta(0,N,1)

plt.plot(msd)
plt.plot(x_fit, y_fit, 'k--')
plt.xlabel("pasos n")
plt.ylabel("<x²>")
plt.title("Desplazamiento cuadrático medio (1D)")
plt.show()


#ahora en dos dimensiones
m = 1000
N = 500

x = np.zeros(m)
y = np.zeros(m)

msd = []

for n in range(N):
    
    pasos_x = np.random.choice([-1,1], size=m)
    pasos_y = np.random.choice([-1,1], size=m)
    
    x += pasos_x
    y += pasos_y
    
    r2 = x**2 + y**2
    msd.append(np.mean(r2))
    
t = np.arange(N)

# vemos cómo encaja con una recta con pendiente 2, es decir, una difusión con D=1/2 ya que <r²>=4Dt
x_fit, y_fit = fit_recta(0,N,2)

plt.plot(msd)
plt.plot(x_fit, y_fit, 'k--')
plt.xlabel("pasos n")
plt.ylabel("<r²>")
plt.title("Desplazamiento cuadrático medio (2D)")
plt.show()

coef = np.polyfit(t, msd, 1)
print("Pendiente =", coef[0])