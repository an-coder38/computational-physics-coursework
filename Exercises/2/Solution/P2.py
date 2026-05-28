import numpy as np
import matplotlib.pyplot as plt

# dv/dt = a - b*v
a=10
b=1
# v should reach the terminal velocity

t_max = 10.0 # time limit
v0 = 0.0 #initial velocity
t_min = 0.0 #starts in 0
dt_list = [0.01, 0.1, 0.5, 1] #time step value

l_dt=len(dt_list)

plt.figure() #opens a figure

for j in range(0, len(dt_list)):
    t = np.arange(t_min, t_max + dt_list[j], dt_list[j]) #equal spaced list of times
    v = np.zeros(len(t)) #matrix of zeros, same length as the time
    v[0] = v0
    for i in range(0,len(t) - 1):
        dv = a-b*v[i] #dv_i but not saved because it is not needed
        v[i+1] = v[i] + dv* dt_list[j]     # Euler algorithm
    plt.plot(t, v, '-', label=rf"Euler $\Delta t$ = {dt_list[j]}")
    
    
# Exact solution solving a EDO: v(t)=(a/b) * (1-exp(-b*t)) (using v0=0)
# from the analytical solution on the problem we observe that the terminal velocity is a/b
t_exact = np.linspace(0, t_max, 200) #creates 200 points of time
v_exact = (a/b)*(1-np.exp(-b*t_exact)) # obtains the exact velocity
plt.plot(t_exact, v_exact, 'k--', label="Solución exacta")

plt.xlabel("t (s)") #label x
plt.ylabel("v (m/s)") #label y
plt.title("Paracaidista: método de Euler") #write a title
plt.legend() #creates a legend
plt.grid() #shows grid
plt.show() #shows the figure