# Código que calcula la trayectoria de proyectiles con distinta 
# celeridad y ángulo inicial, teniendo en cuenta la resistencia del aire.
import numpy as np
import matplotlib.pyplot as plt
from math import cos, sin, sqrt

# aprox adiab solo válida cuando lo del paréntesis es >0
# la variabilidad de la gravedad con la altura no afecta considerablemente a ningún modelo
# a alturas bajas pero, salvo en el primer modelo, se considerará variable para no perder generalidad

#para obtener todas las gráficas de altura en función de alcance para cada modelo se puede descomentar la parte de visualización de la función modelo

# constantes que van a usar los modelos
g = 9.81 # aceleración de la gravedad
R_T = 6300*1000 # m, radio de la tierra
B2m = 4e-5 #constantes
alt_0 = 10**4 # m
a = 6.5e-3 # K/m
alpha = 2.5
T_0 = 280 # K, temperatura a nivel del mar


def modelo(t0, dt, x0, y0, v0, theta_in, theta_fin, dtheta, modo):
    """
    modo:
    0 -> En vacío con g cte
    1 -> En vacío con g variable
    2 -> Aire con densidad constante
    3 -> Aproximación isotérmica
    4 -> Aproximación adiabática
    """

    nombres = [
        'Modelo sin aire con g cte',
        'Modelo sin aire con g variable',
        'Modelo con densidad de aire cte',
        'Modelo con aproximación isotérmica',
        'Modelo con aproximación adiabática'
    ]

    print(nombres[modo])

    theta_degree = np.arange(theta_in, theta_fin + dtheta, dtheta)
    theta = np.radians(theta_degree)

    alcance = []
#    plt.figure()

    for j in range(len(theta)):
        # Bucle sobre todos los ángulos
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

            # g constante o variable
            if modo == 0:
                g_eff = g
            else:
                g_eff = g * (R_T / (R_T + y[i]))**2

            # términos de rozamiento
            if modo == 0 or modo == 1:
                drag_x = 0
                drag_y = 0

            elif modo == 2:
                drag_x = B2m * v[i] * vx[i]
                drag_y = B2m * v[i] * vy[i]

            elif modo == 3:
                factor = np.exp(-y[i] / alt_0)
                drag_x = B2m * v[i] * vx[i] * factor
                drag_y = B2m * v[i] * vy[i] * factor

            elif modo == 4:
                factor = (1 - a * y[i] / T_0)**alpha
                drag_x = B2m * v[i] * vx[i] * factor
                drag_y = B2m * v[i] * vy[i] * factor

            vx.append(vx[i] - drag_x * dt)
            vy.append(vy[i] - g_eff * dt - drag_y * dt)
            v.append(sqrt(vx[i]**2 + vy[i]**2))
            i += 1

        # alcance aproximado válido para dt pequeños
        alcance.append(x[-2] + (x[-1] - x[-2]) / 2)

#        plt.plot(x, y, linewidth=1, label=f"theta={np.degrees(theta[j]):.1f}º")
#        print(f"{np.degrees(theta[j]):.1f}  {alcance[j]:.1f}")

    # cálculo del alcance máximo
    alcance_maximo = alcance[0]
    i_max = 0
    for i in range(len(alcance)):
        if alcance[i] > alcance_maximo:
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
t0 = 0; dt = 0.1
x0 = 0; y0 = 0
v0 = 700
dtheta = 1; theta_in =0; theta_fin = 90 #condiciones iniciales del barrido angular
theta_degree = np.arange(theta_in, theta_fin+dtheta, dtheta)

# calculo x, y, vx, vy, v, alcance y el alcande máximo junto a su índice correspondiente para luego poder seleccionar el ángulo que nos da el máximo alcance
x_0, y_0, vx_0, vy_0, v_0, alcance_0, alcance_maximo_0, i0 = modelo(t0, dt, x0, y0, v0, theta_in, theta_fin, dtheta, 0)
print(f'Alcance máximo: {alcance_maximo_0:.1f} m con ángulo: {theta_degree[i0]:.1f}º')
x_g, y_g, vx_g, vy_g, v_g, alcance_g, alcance_maximo_g, ig = modelo(t0, dt, x0, y0, v0, theta_in, theta_fin, dtheta, 1)
print(f'Alcance máximo: {alcance_maximo_g:.1f} m con ángulo: {theta_degree[ig]:.1f}º')
x_a, y_a, vx_a, vy_a, v_a, alcance_a, alcance_maximo_a, ia = modelo(t0, dt, x0, y0, v0, theta_in, theta_fin, dtheta, 2)
print(f'Alcance máximo: {alcance_maximo_a:.1f} m con ángulo: {theta_degree[ia]:.1f}º')
x_t, y_t, vx_t, vy_t, v_t, alcance_t, alcance_maximo_t, it = modelo(t0, dt, x0, y0, v0, theta_in, theta_fin, dtheta, 3)
print(f'Alcance máximo: {alcance_maximo_t:.1f} m con ángulo: {theta_degree[it]:.1f}º')
x_ad, y_ad, vx_ad, vy_ad, v_ad, alcance_ad, alcance_maximo_ad, iad = modelo(t0, dt, x0, y0, v0, theta_in, theta_fin, dtheta, 4)
print(f'Alcance máximo: {alcance_maximo_ad:.1f} m con ángulo: {theta_degree[iad]:.1f}º')

# creo figura para comparar los alcances de los diferentes métodos
plt.figure()
plt.title('Comparación de alcances de los diferentes métodos')
plt.plot(theta_degree, alcance_0, label='En vacío, g cte')
plt.plot(theta_degree, alcance_g, label='En vacío, g variable')
plt.plot(theta_degree, alcance_a, label='Con aire, densidad cte')
plt.plot(theta_degree, alcance_t, label='Aprox isotérmica')
plt.plot(theta_degree, alcance_ad, label='Aprox adiabática')
plt.legend()
plt.xlabel("theta(º)")             # Nombre a los ejes
plt.ylabel("Alcance (m)")
plt.title('Comparación de modelos')
plt.show()  
plt.close()
#print(f"El alcance máximo es de {alcance_maximo:.1f} y se alcanza a los {theta_degree[i_max]:.1f} grados")


#Creo ahora el caso de alcance máximo para cada modelo, tomando solo el caso para el theta que nos ofrece el mayor alcance en cada modelo
x_0_max, y_0_max, *_ = modelo(t0, dt, x0, y0, v0, theta_degree[i0], theta_degree[i0], dtheta, 0)
x_g_max, y_g_max, *_ = modelo(t0, dt, x0, y0, v0, theta_degree[ig], theta_degree[ig], dtheta, 1)
x_a_max, y_a_max, *_ = modelo(t0, dt, x0, y0, v0, theta_degree[ia], theta_degree[ia], dtheta, 2)
x_t_max, y_t_max, *_ = modelo(t0, dt, x0, y0, v0, theta_degree[it], theta_degree[it], dtheta, 3)
x_ad_max, y_ad_max, *_ = modelo(t0, dt, x0, y0, v0, theta_degree[iad], theta_degree[iad], dtheta, 4)

# creo figura para comparar los alcances máximos de los diferentes métodos para las mismas condiciones iniciales
plt.figure()
plt.title('Comparación de alcances de los diferentes métodos')
plt.plot(x_0_max, y_0_max, label='En vacío, g cte')
plt.plot(x_g_max, y_g_max, label='En vacío, g variable')
plt.plot(x_a_max, y_a_max, label='Con aire, densidad cte')
plt.plot(x_t_max, y_t_max, label='Aprox isotérmica')
plt.plot(x_ad_max, y_ad_max, label='Aprox adiabática')
plt.legend()
plt.xlabel("Alcance (m)")             # Nombre a los ejes
plt.ylabel("Altura (m)")
plt.title('Comparación de modelos')
plt.show()  
plt.close()
