import numpy as np
import matplotlib.pyplot as plt

#initial parameters
tau_a = 1
tau_b = 1.5
na0 = 25
nb0 = 5

t_max = 10.0 # time limit
t_min = 0.0 #starts in 0
dt_list = [0.01] #time step value

l_dt=len(dt_list)

plt.figure() #opens a figure

for j in range(0, len(dt_list)):
    t = np.arange(t_min, t_max + dt_list[j], dt_list[j]) #equal spaced list of times
    na = np.zeros(len(t)) #matrix of zeros, same length as the time
    nb = np.zeros(len(t))
    na[0] = na0
    nb[0] = nb0
    for i in range(0,len(t) - 1):
        dna = -na[i]/tau_a
        na[i+1] = na[i] + dna*dt_list[j]     # Euler algorithm
        dnb = na[i]/tau_a -nb[i]/tau_b
        nb[i+1] = nb[i] + dnb*dt_list[j]
    plt.plot(t, na, '-', label=rf"Na -- Euler $\Delta t$ = {dt_list[j]}")
    plt.plot(t, nb, '-', label=rf"Nb -- Euler $\Delta t$ = {dt_list[j]}")
    

# Exact solution solving both EDOs analytically
t_exact = np.linspace(0, t_max, 200) #creates 200 points of time
na_exact = na0*(np.exp(-t_exact/tau_a)) # obtains the exact na
plt.plot(t_exact, na_exact, 'k--', label="Exacta na")

k=tau_b*na0/(tau_a-tau_b)
nb_exact = (nb0-k)*(np.exp(-t_exact/tau_b))+k*(np.exp(-t_exact/tau_a)) # obtains the exact nb
plt.plot(t_exact, nb_exact, 'k--', label="Exacta nb")


plt.xlabel("t (s)") #label x
plt.ylabel("n (áts)") #label y
plt.title("Decaimiento de núcleos atómicos: Núcleo padre y primer hijo") #write a title
plt.legend() #creates a legend
plt.grid() #shows grid
plt.show() #shows the figure

# tau_b much bigger than tau_a, most of the nucleus turn into b and after some time it decays into c (another nucleus)
#tau_a much bigger than tau_b, na(t)~na0 and it exponencially decreases turning into b that quickly turns into c (a decreases exponentilly but no b rapidly decays into c)
# tau_a near tau_b, both decays compete and it ends up decaying after some time
# tau_a=tau_b, error in the division by 0 so the analytical cannot be painted