import numpy as np
import matplotlib.pyplot as plt
from math import cos, sin, sqrt

# constantes que van a usar los modelos
g = 9.81 # aceleración de la gravedad
R_T = 6300*1000 # m, radio de la tierra
B2m = 4e-5 #constantes
alt_0 = 10**4 # m
a = 6.5e-3 # K/m
alpha = 2.5
T_0 = 280 # K, temperatura a nivel del mar
vd = 35 #m/s
delta = 5 #m/s


def modelo(t0, dt, x0, y0, v0, theta_in, theta_fin, dtheta, vwx, vwy, modo):
    """
    modo:
    0 -> Sin viento # caso equivalente a vwx=vwy=0
    1 -> Con viento 
    """
    nombres = ['Modelo sin viento', 'Modelo con viento']
    print(nombres[modo])
    theta_degree = np.arange(theta_in, theta_fin + dtheta, dtheta)
    theta = np.radians(theta_degree)

    alcance = []
#    plt.figure()

    for j in range(len(theta)): # Bucle sobre todos los ángulos
        
        t, x, y, vx, vy, v = [], [], [], [], [], []

        # Condiciones iniciales
        t.append(t0)
        x.append(x0)
        y.append(y0)
        vx.append(v0 * cos(theta[j]))
        vy.append(v0 * sin(theta[j]))
        v.append(v0)

        i = 0
        while y[i] >= 0:
            t.append(t[i] + dt)
            x.append(x[i] + vx[i] * dt)
            y.append(y[i] + vy[i] * dt)

            g_eff = g * (R_T / (R_T + y[i]))**2 #aunque de forma poco relevante en alturas bajas, la gravedad depende la altura

            # término de rozamiento
            B2m = 0.0039 + 0.0058/(1+np.exp((v[i]-vd)/delta))

            if modo == 0:
                drag_x = B2m * v[i] * vx[i]
                drag_y = B2m * v[i] * vy[i]

            elif modo == 1:
                # velocidad relativa de la bola respecto al viento
                v_con_x = vx[i] - vwx 
                v_con_y = vy[i] - vwy 
                v_con = np.sqrt((v_con_x)**2 + (v_con_y)**2) # |v+vw| con v, vw vectores
                
                drag_x = B2m * v_con * v_con_x
                drag_y = B2m * v_con * v_con_y

            vx.append(vx[i] - drag_x * dt)
            vy.append(vy[i] - g_eff * dt - drag_y * dt)
            v.append(sqrt(vx[i+1]**2 + vy[i+1]**2))
            i += 1

        # alcance aproximado válido para dt pequeños
        alcance.append((x[-1] + x[-2]) / 2)

#        plt.plot(x, y, linewidth=1, label=f"theta={np.degrees(theta[j]):.1f}º")
#        print(f"{np.degrees(theta[j]):.1f}  {alcance[j]:.1f}")

    # identificación del alcance máximo, precisión dtheta
    alcance_maximo = alcance[0]
    i_max = 0
    for i in range(len(alcance)):
        if alcance[i] > alcance_maximo: # busca el índice del alcance máximo
            alcance_maximo = alcance[i]
            i_max = i

#    plt.xlabel("Alcance (m)")
#    plt.ylabel("Altura (m)")
#    plt.title(nombres[modo])
    # plt.legend()  # evitar usarla si dtheta es bajo
#    plt.show()
#    plt.close()

#    plt.figure()
#    plt.plot(theta_degree, alcance)
#    plt.xlabel("theta(º)")
#    plt.ylabel("Alcance (m)")
#    plt.title(nombres[modo])
#    plt.show()
#    plt.close()

    return x, y, vx, vy, v, alcance, alcance_maximo, i_max



# Condiciones iniciales del problema
t0 = 0; dt = 0.01
x0 = 0; y0 = 0
v0 = 49
vwx = 49/11; vwy = 0 # velocidad del viento
dtheta = 0.1; theta_in =0; theta_fin = 90 #condiciones iniciales del barrido angular
theta_degree = np.arange(theta_in, theta_fin+dtheta, dtheta)

x_0, y_0, vx_0, vy_0, v_0, alcance_0, alcance_maximo_0, i0 = modelo(t0, dt, x0, y0, v0, theta_in, theta_fin, dtheta, vwx, vwy, 0)
print(f'Alcance máximo: {alcance_maximo_0:.1f} m con ángulo: {theta_degree[i0]:.1f}º')
x_head, y_head, vx_head, vy_head, v_head, alcance_head, alcance_maximo_head, ihead = modelo(t0, dt, x0, y0, v0, theta_in, theta_fin, dtheta, vwx, vwy, 1)
print(f'Alcance máximo con viento a favor: {alcance_maximo_head:.1f} m con ángulo: {theta_degree[ihead]:.1f}º')
x_tail, y_tail, vx_tail, vy_tail, v_tail, alcance_tail, alcance_maximo_tail, itail = modelo(t0, dt, x0, y0, v0, theta_in, theta_fin, dtheta, -vwx, -vwy, 1)
print(f'Alcance máximo con viento en contra: {alcance_maximo_tail:.1f} m con ángulo: {theta_degree[itail]:.1f}º')

# creo figura para comparar los alcances de los diferentes métodos
plt.figure()
plt.plot(theta_degree, alcance_0, label='Sin viento')
plt.plot(theta_degree, alcance_head, label='Con viento a favor')
plt.plot(theta_degree, alcance_tail, label='Con viento en contra')
plt.legend()
plt.xlabel("theta(º)")             # Nombre a los ejes
plt.ylabel("Alcance máx (m)")
# plt.title('Alcance de una pelota de béisbol en función del ángulo')
plt.savefig("beisbol_alcance.png", dpi=300, bbox_inches="tight")
plt.show()  
plt.close()
#print(f"El alcance máximo es de {alcance_maximo:.1f} y se alcanza a los {theta_degree[i_max]:.1f} grados")


#Creo ahora el caso de alcance máximo para cada modelo, tomando solo el caso para el theta que nos ofrece el mayor alcance en cada modelo
x_0_max, y_0_max, *_ = modelo(t0, dt, x0, y0, v0, theta_degree[i0], theta_degree[i0], dtheta, vwx, vwy, 0)
x_head_max, y_head_max, *_ = modelo(t0, dt, x0, y0, v0, theta_degree[ihead], theta_degree[ihead], dtheta, vwx, vwy, 1)
x_tail_max, y_tail_max, *_ = modelo(t0, dt, x0, y0, v0, theta_degree[itail], theta_degree[itail], dtheta, -vwx, -vwy, 1)

# creo figura para comparar los alcances máximos de los diferentes métodos para las mismas condiciones iniciales
plt.figure()
plt.plot(x_0_max, y_0_max, label='Sin viento')
plt.plot(x_head_max, y_head_max, label='Con viento a favor')
plt.plot(x_tail_max, y_tail_max, label='Con viento en contra')
# plt.legend()
plt.xlabel("Alcance (m)")             # Nombre a los ejes
plt.ylabel("Altura (m)")
# plt.title('Lanzamiento optimizado de una pelota de béisbol')
plt.savefig("beisbol_optimizado.png", dpi=300, bbox_inches="tight")
plt.show()  
plt.close()

# Creo ahora el caso de alcance máximo para cada modelo, tomando solo un ángulo en concreto
theta=89 # en este caso lo dejo en 89º para observar cómo el viento en contra puede ser altamente relevante e incluso invertir la dirección de la trayectoria
x_0_max, y_0_max, *_ = modelo(t0, dt, x0, y0, v0, theta, theta, dtheta, vwx, vwy, 0)
x_head_max, y_head_max, *_ = modelo(t0, dt, x0, y0, v0, theta, theta, dtheta, vwx, vwy, 1)
x_tail_max, y_tail_max, *_ = modelo(t0, dt, x0, y0, v0, theta, theta, dtheta, -vwx, -vwy, 1)

# creo figura para comparar las trayectorias de los diferentes métodos para un mismo ángulo
plt.figure()
plt.plot(x_0_max, y_0_max, label='Sin viento')
plt.plot(x_head_max, y_head_max, label='Con viento a favor')
plt.plot(x_tail_max, y_tail_max, label='Con viento en contra')
# plt.legend()
plt.xlabel("Alcance (m)")             # Nombre a los ejes
plt.ylabel("Altura (m)")
# plt.title(f'Lanzamiento pelota de béisbol: ángulo inicial: {theta:.0f}º')
plt.savefig("beisbol_trayectorias.png", dpi=300, bbox_inches="tight")
plt.show()  
plt.close()


