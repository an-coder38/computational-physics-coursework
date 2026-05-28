import numpy as np
import matplotlib.pyplot as plt
from math import pi
# Constantes físicas del problema
g = 9.81 # aceleración gravitatoria
rho = 1.2 # kg/m^3, densidad del aire
m = 45.93/1000 # kg, masa pelota de golf
R = 0.02135 # m, radio pelota de golf
A = pi*R**2 # área pelota de golf
R_T = 6300*1000 # m, radio de la tierra

# Programa que utilizando las funciones definidas, calcula la trayectoria bajo 
# varias condiciones iniciales tanto con viento como sin viento y luego obtiene 
# las gráficas más esclarecedoras para analizar la trayectorias de los lanzamientos

# Modelo 3D con Magnus
def modelo_golf_3D(dt, r, v0, theta_deg, phi_deg, mSomega_vec, vw, rugosidad):
    
    # añado un parámetro rugosidad: 0 si lisa, 1 si rugosa
    # para el rango de velocidades del problema de golf se puede considerar cte
    # el valor de C=1 para la fuerza de arrastre si es lisa
    
    theta = np.radians(theta_deg)   # ángulo vertical
    phi = np.radians(phi_deg)       # ángulo horizontal

    # velocidad inicial en 3D desde coord esféricas
    vx = v0 * np.cos(theta) * np.cos(phi)
    vy = v0 * np.cos(theta) * np.sin(phi)
    vz = v0 * np.sin(theta)
    
    v = np.array([vx, vy, vz])
    
    trayectoria = [r.copy()]
    
    while r[2] >= 0: #z>=0
        
        v_rel = v - vw # tengo en cuenta el efecto del viento
        v_mod = np.linalg.norm(v_rel) # módulo de la velocidad
        
        # C=1 si v<14; C=14/v si v>14
        #if v_mod < 14:
        #    C=1
        #else:
        #    C=14/v_mod
        if rugosidad == 0: C=1 # si lisa
        if rugosidad == 1: C = 14/max(v_mod, 14) # si rugosa, evita un bucle if extra
        
        # suponemos densidad del aire aproximadamente constante
        a_drag = -0.5* C * rho* A * v_mod * v_rel /m # contribución debida al drag
        
        a_magnus = np.cross(mSomega_vec, v_rel) # Magnus, obtenido con un producto vectorial
        
        g_eff = g * (R_T / (R_T + r[2]))**2  # la gravedad depende la altura
        a_grav = np.array([0, 0, -g_eff]) # contribución gravitatoria, que no sea cte es poco relevante en este problema
        
        a = a_grav + a_drag + a_magnus # aceleración total
        
        # Euler
        r = r + v * dt 
        v = v + a * dt
        trayectoria.append(r.copy())
    
    return np.array(trayectoria)


def optimizar_theta(dt, r0, v0, theta_in, theta_fin, dtheta, phi_deg, msomega_vec, vw, rugosidad):
    
    theta_degree = np.arange(theta_in, theta_fin + dtheta, dtheta)
    alcances = []

    for theta in theta_degree:
        tray = modelo_golf_3D(dt, r0, v0, theta, phi_deg, msomega_vec, vw, rugosidad)
        alcance = (tray[-2,0] + tray [-1,0])/2
        alcances.append(alcance)

    #alcances = np.array(alcances)
    i_max = np.argmax(alcances)

    return theta_degree[i_max], alcances[i_max]



# Condiciones iniciales del problema
dt = 0.01
v0 = 70 # velocidad inicial de la pelota de golf
theta = 12 # ángulo theta cualquiera, luego se optimizará
phi = 0 # fijo la dirección
r0 = np.array([0.0, 0.0, 0.0]) # posición inicial
msomega = 0.25  # rad/s
msomega_hook = np.array([0, 0, -msomega]) # hook en dir -z
msomega_slice = np.array([0, 0, msomega]) # slice en dir z

vw = np.array([3,4,5]) # velocidad del viento

#Comenzamos con un caso no optimizado

# PELOTA RUGOSA

tray_hook = modelo_golf_3D(dt, r0, v0, theta, phi, msomega_hook, vw, 1) #cálculo de trayectorias
tray_slice = modelo_golf_3D(dt, r0, v0, theta, phi, msomega_slice, vw, 1)
tray_no = modelo_golf_3D(dt, r0, v0, theta, phi, 0*msomega_slice, vw, 1)

#Mostrar en la terminal resultados interesantes para el estudio de los lanzamientos
print('Pelota de golf rugosa')
alcance_hook = (tray_hook[-2,0]+tray_hook[-1,0])/2
print(f"Desviación lateral hook: {abs(tray_hook[-2,1]+tray_hook[-1,1])/2:.1f}")
print(f"Alcance hook: {alcance_hook:.1f}")

alcance_slice = (tray_slice[-2,0]+tray_slice[-1,0])/2
print(f"Desviación lateral slice: {abs(tray_slice[-2,1]+tray_slice[-1,1])/2:.1f}")
print(f"Alcance slice: {alcance_slice:.1f}")

alcance_no = (tray_no[-2,0]+tray_no[-1,0])/2
print(f"Desviación lateral sin Magnus: {abs(tray_no[-2,1]+tray_no[-1,1])/2:.1f}")
print(f"Alcance sin Magnus: {alcance_no:.1f}")

# PELOTA LISA: cálculo de trayectorias y, alcance y desviación lateral

tray_hookl = modelo_golf_3D(dt, r0, v0, theta, phi, msomega_hook, vw, 0)
tray_slicel = modelo_golf_3D(dt, r0, v0, theta, phi, msomega_slice, vw, 0)
tray_nol = modelo_golf_3D(dt, r0, v0, theta, phi, 0*msomega_slice, vw, 0)

print('Pelota de golf lisa')
alcance_hookl = (tray_hookl[-2,0]+tray_hookl[-1,0])/2
print(f"Desviación lateral hook: {abs(tray_hookl[-2,1]+tray_hookl[-1,1])/2:.1f}")
print(f"Alcance hook: {alcance_hookl:.1f}")

alcance_slicel = (tray_slicel[-2,0]+tray_slicel[-1,0])/2
print(f"Desviación lateral slice: {abs(tray_slicel[-2,1]+tray_slicel[-1,1])/2:.1f}")
print(f"Alcance slice: {alcance_slicel:.1f}")

alcance_nol = (tray_nol[-2,0]+tray_nol[-1,0])/2
print(f"Desviación lateral sin Magnus: {abs(tray_nol[-2,1]+tray_nol[-1,1])/2:.1f}")
print(f"Alcance sin Magnus: {alcance_nol:.1f}")

# Optimizar el ángulo y repetir el procedimiento

theta_in=0; theta_fin=90; dtheta=1;

print(f"\nOPTIMIZACIÓN DEL ÁNGULO (precisión {dtheta}°)\n")

# --- RUGOSA ---

theta_opt_hook_r, alcance_opt_hook_r = optimizar_theta(dt, r0, v0, theta_in, theta_fin, dtheta, phi, msomega_hook, vw, 1)
theta_opt_slice_r, alcance_opt_slice_r = optimizar_theta(dt, r0, v0, theta_in, theta_fin, dtheta, phi, msomega_slice, vw, 1)
theta_opt_no_r, alcance_opt_no_r = optimizar_theta(dt, r0, v0, theta_in, theta_fin, dtheta, phi, 0*msomega_slice, vw, 1)

# --- LISA ---

theta_opt_hook_l, alcance_opt_hook_l = optimizar_theta(dt, r0, v0, theta_in, theta_fin, dtheta, phi, msomega_hook, vw, 0)
theta_opt_slice_l, alcance_opt_slice_l = optimizar_theta(dt, r0, v0, theta_in, theta_fin, dtheta, phi, msomega_slice, vw, 0)
theta_opt_no_l, alcance_opt_no_l = optimizar_theta(dt, r0, v0, theta_in, theta_fin, dtheta, phi, 0*msomega_slice, vw, 0)

# Sobreescribo los resultados de las trayectorias obtenidas para un ángulo cualquiera con los datos de las trayectorias que maximizan el alcance en función de theta inicial

# PELOTA RUGOSA optimizada

print("Pelota rugosa")
print(f"Hook → theta óptimo = {theta_opt_hook_r}°, alcance = {alcance_opt_hook_r:.1f} m")
print(f"Slice → theta óptimo = {theta_opt_slice_r}°, alcance = {alcance_opt_slice_r:.1f} m")
print(f"Sin Magnus → theta óptimo = {theta_opt_no_r}°, alcance = {alcance_opt_no_r:.1f} m")

tray_hook = modelo_golf_3D(dt, r0, v0, theta_opt_hook_r, phi, msomega_hook, vw, 1)
tray_slice = modelo_golf_3D(dt, r0, v0, theta_opt_slice_r, phi, msomega_slice, vw, 1)
tray_no = modelo_golf_3D(dt, r0, v0, theta_opt_no_r, phi, 0*msomega_slice, vw, 1)

print('Pelota de golf rugosa')
alcance_hook = (tray_hook[-2,0]+tray_hook[-1,0])/2
print(f"Desviación lateral hook: {abs(tray_hook[-2,1]+tray_hook[-1,1])/2:.1f}")
print(f"Alcance hook: {alcance_hook:.1f}")

alcance_slice = (tray_slice[-2,0]+tray_slice[-1,0])/2
print(f"Desviación lateral slice: {abs(tray_slice[-2,1]+tray_slice[-1,1])/2:.1f}")
print(f"Alcance slice: {alcance_slice:.1f}")

alcance_no = (tray_no[-2,0]+tray_no[-1,0])/2
print(f"Desviación lateral sin Magnus: {abs(tray_no[-2,1]+tray_no[-1,1])/2:.1f}")
print(f"Alcance sin Magnus: {alcance_no:.1f}")

# PELOTA LISA optimizada

print("\nPelota lisa")
print(f"Hook → theta óptimo = {theta_opt_hook_l}°, alcance = {alcance_opt_hook_l:.1f} m")
print(f"Slice → theta óptimo = {theta_opt_slice_l}°, alcance = {alcance_opt_slice_l:.1f} m")
print(f"Sin Magnus → theta óptimo = {theta_opt_no_l}°, alcance = {alcance_opt_no_l:.1f} m")

tray_hookl = modelo_golf_3D(dt, r0, v0, theta_opt_hook_l, phi, msomega_hook, vw, 0)
tray_slicel = modelo_golf_3D(dt, r0, v0, theta_opt_slice_l, phi, msomega_slice, vw, 0)
tray_nol = modelo_golf_3D(dt, r0, v0, theta_opt_no_l, phi, 0*msomega_slice, vw, 0)

print('Pelota de golf lisa')
alcance_hookl = (tray_hookl[-2,0]+tray_hookl[-1,0])/2
print(f"Desviación lateral hook: {abs(tray_hookl[-2,1]+tray_hookl[-1,1])/2:.1f}")
print(f"Alcance hook: {alcance_hookl:.1f}")

alcance_slicel = (tray_slicel[-2,0]+tray_slicel[-1,0])/2
print(f"Desviación lateral slice: {abs(tray_slicel[-2,1]+tray_slicel[-1,1])/2:.1f}")
print(f"Alcance slice: {alcance_slicel:.1f}")

alcance_nol = (tray_nol[-2,0]+tray_nol[-1,0])/2
print(f"Desviación lateral sin Magnus: {abs(tray_nol[-2,1]+tray_nol[-1,1])/2:.1f}")
print(f"Alcance sin Magnus: {alcance_nol:.1f}")

# Creación de figuras para los ángulos optimizados

# Comparación de lisa y rugosa para hook y slice
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.plot(tray_hook[:,0], tray_hook[:,1], tray_hook[:,2], label = 'hook -- rugosa')
ax.plot(tray_slice[:,0], tray_slice[:,1], tray_slice[:,2], label = 'slice -- rugosa')
ax.plot(tray_no[:,0], tray_no[:,1], tray_no[:,2], label = 'no Magnus -- rugosa')
ax.plot(tray_hookl[:,0], tray_hookl[:,1], tray_hookl[:,2], label = 'hook -- lisa')
ax.plot(tray_slicel[:,0], tray_slicel[:,1], tray_slicel[:,2], label = 'slice -- lisa')
ax.plot(tray_nol[:,0], tray_nol[:,1], tray_nol[:,2], label = 'no Magnus -- lisa')
# plt.legend()
ax.set_xlabel("X (m)")
ax.set_ylabel("Y (m)")
ax.set_zlabel("Z (m)")
# ax.set_title("Pelota de golf rugosa vs lisa")
plt.savefig("golf_vscon.png", dpi=300)
plt.show()

# Vista desde un punto de vista superior del plano xy
plt.figure()
plt.plot(tray_hook[:,0], tray_hook[:,1], label='hook -- rugosa')
plt.plot(tray_slice[:,0], tray_slice[:,1], label='slice -- rugosa')
plt.plot(tray_no[:,0], tray_no[:,1], label='no Magnus -- rugosa')
plt.plot(tray_hookl[:,0], tray_hookl[:,1], label='hook -- lisa')
plt.plot(tray_slicel[:,0], tray_slicel[:,1], label='slice -- lisa')
plt.plot(tray_nol[:,0], tray_nol[:,1], label='no Magnus -- lisa')
plt.xlabel("X (m)")
plt.ylabel("Y (m)")
# plt.legend()
plt.savefig("golf_lateralXYcon.png", dpi=300, bbox_inches="tight")
# plt.title("Vista superior (desviación lateral)")
plt.show()

# Vista desde un punto de vista superior del plano xz
plt.figure()
plt.plot(tray_hook[:,0], tray_hook[:,2], label='hook -- rugosa')
plt.plot(tray_slice[:,0], tray_slice[:,2], label='slice -- rugosa')
plt.plot(tray_no[:,0], tray_no[:,2], label='no Magnus -- rugosa')
plt.plot(tray_hookl[:,0], tray_hookl[:,2], label='hook -- lisa')
plt.plot(tray_slicel[:,0], tray_slicel[:,2], label='slice -- lisa')
plt.plot(tray_nol[:,0], tray_nol[:,2], label='no Magnus -- lisa')
plt.xlabel("X (m)")
plt.ylabel("Z (m)")
# plt.legend()
# plt.title("Vista lateral (XZ)")
plt.savefig("golf_lateralXZcon.png", dpi=300, bbox_inches="tight")
plt.show()

# Vista desde un punto de vista superior del plano yz
plt.figure()
plt.plot(tray_hook[:,1], tray_hook[:,2], label='hook -- rugosa')
plt.plot(tray_slice[:,1], tray_slice[:,2], label='slice -- rugosa')
plt.plot(tray_no[:,1], tray_no[:,2], label='no Magnus -- rugosa')
plt.plot(tray_hookl[:,1], tray_hookl[:,2], label='hook -- lisa')
plt.plot(tray_slicel[:,1], tray_slicel[:,2], label='slice -- lisa')
plt.plot(tray_nol[:,1], tray_nol[:,2], label='no Magnus -- lisa')
plt.xlabel("Y (m)")
plt.ylabel("Z (m)")
# plt.legend()
# plt.title("Vista lateral (YZ)")
plt.savefig("golf_lateralYZcon.png", dpi=300, bbox_inches="tight")
plt.show()

###################CASO SIN VIENTO#######################################

vw = np.array([0,0,0]) # velocidad del viento

print('CASO SIN VIENTO')

# Optimizar el ángulo

theta_in=0; theta_fin=90; dtheta=1;

print(f"\nOPTIMIZACIÓN DEL ÁNGULO (precisión {dtheta}°)\n")

# --- RUGOSA ---

theta_opt_hook_r, alcance_opt_hook_r = optimizar_theta(dt, r0, v0, theta_in, theta_fin, dtheta, phi, msomega_hook, vw, 1)
theta_opt_slice_r, alcance_opt_slice_r = optimizar_theta(dt, r0, v0, theta_in, theta_fin, dtheta, phi, msomega_slice, vw, 1)
theta_opt_no_r, alcance_opt_no_r = optimizar_theta(dt, r0, v0, theta_in, theta_fin, dtheta, phi, 0*msomega_slice, vw, 1)

# --- LISA ---

theta_opt_hook_l, alcance_opt_hook_l = optimizar_theta(dt, r0, v0, theta_in, theta_fin, dtheta, phi, msomega_hook, vw, 0)
theta_opt_slice_l, alcance_opt_slice_l = optimizar_theta(dt, r0, v0, theta_in, theta_fin, dtheta, phi, msomega_slice, vw, 0)
theta_opt_no_l, alcance_opt_no_l = optimizar_theta(dt, r0, v0, theta_in, theta_fin, dtheta, phi, 0*msomega_slice, vw, 0)

# Sobreescribo los resultados de las trayectorias obtenidas para un ángulo cualquiera con los datos de las trayectorias que maximizan el alcance en función de theta inicial

# PELOTA RUGOSA optimizada

print("Pelota rugosa")
print(f"Hook → theta óptimo = {theta_opt_hook_r}°, alcance = {alcance_opt_hook_r:.1f} m")
print(f"Slice → theta óptimo = {theta_opt_slice_r}°, alcance = {alcance_opt_slice_r:.1f} m")
print(f"Sin Magnus → theta óptimo = {theta_opt_no_r}°, alcance = {alcance_opt_no_r:.1f} m")

tray_hook = modelo_golf_3D(dt, r0, v0, theta_opt_hook_r, phi, msomega_hook, vw, 1)
tray_slice = modelo_golf_3D(dt, r0, v0, theta_opt_slice_r, phi, msomega_slice, vw, 1)
tray_no = modelo_golf_3D(dt, r0, v0, theta_opt_no_r, phi, 0*msomega_slice, vw, 1)

print('Pelota de golf rugosa')
alcance_hook = (tray_hook[-2,0]+tray_hook[-1,0])/2
print(f"Desviación lateral hook: {abs(tray_hook[-2,1]+tray_hook[-1,1])/2:.1f}")
print(f"Alcance hook: {alcance_hook:.1f}")

alcance_slice = (tray_slice[-2,0]+tray_slice[-1,0])/2
print(f"Desviación lateral slice: {abs(tray_slice[-2,1]+tray_slice[-1,1])/2:.1f}")
print(f"Alcance slice: {alcance_slice:.1f}")

alcance_no = (tray_no[-2,0]+tray_no[-1,0])/2
print(f"Desviación lateral sin Magnus: {abs(tray_no[-2,1]+tray_no[-1,1])/2:.1f}")
print(f"Alcance sin Magnus: {alcance_no:.1f}")

# PELOTA LISA optimizada

print("\nPelota lisa")
print(f"Hook → theta óptimo = {theta_opt_hook_l}°, alcance = {alcance_opt_hook_l:.1f} m")
print(f"Slice → theta óptimo = {theta_opt_slice_l}°, alcance = {alcance_opt_slice_l:.1f} m")
print(f"Sin Magnus → theta óptimo = {theta_opt_no_l}°, alcance = {alcance_opt_no_l:.1f} m")

tray_hookl = modelo_golf_3D(dt, r0, v0, theta_opt_hook_l, phi, msomega_hook, vw, 0)
tray_slicel = modelo_golf_3D(dt, r0, v0, theta_opt_slice_l, phi, msomega_slice, vw, 0)
tray_nol = modelo_golf_3D(dt, r0, v0, theta_opt_no_l, phi, 0*msomega_slice, vw, 0)

print('Pelota de golf lisa')
alcance_hookl = (tray_hookl[-2,0]+tray_hookl[-1,0])/2
print(f"Desviación lateral hook: {abs(tray_hookl[-2,1]+tray_hookl[-1,1])/2:.1f}")
print(f"Alcance hook: {alcance_hookl:.1f}")

alcance_slicel = (tray_slicel[-2,0]+tray_slicel[-1,0])/2
print(f"Desviación lateral slice: {abs(tray_slicel[-2,1]+tray_slicel[-1,1])/2:.1f}")
print(f"Alcance slice: {alcance_slicel:.1f}")

alcance_nol = (tray_nol[-2,0]+tray_nol[-1,0])/2
print(f"Desviación lateral sin Magnus: {abs(tray_nol[-2,1]+tray_nol[-1,1])/2:.1f}")
print(f"Alcance sin Magnus: {alcance_nol:.1f}")

# Comparación de lisa y rugosa para hook y slice
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.plot(tray_hook[:,0], tray_hook[:,1], tray_hook[:,2], label = 'hook -- rugosa')
ax.plot(tray_slice[:,0], tray_slice[:,1], tray_slice[:,2], label = 'slice -- rugosa')
ax.plot(tray_no[:,0], tray_no[:,1], tray_no[:,2], label = 'no Magnus -- rugosa')
ax.plot(tray_hookl[:,0], tray_hookl[:,1], tray_hookl[:,2], label = 'hook -- lisa')
ax.plot(tray_slicel[:,0], tray_slicel[:,1], tray_slicel[:,2], label = 'slice -- lisa')
ax.plot(tray_nol[:,0], tray_nol[:,1], tray_nol[:,2], label = 'no Magnus -- lisa')
# plt.legend()
ax.set_xlabel("X (m)")
ax.set_ylabel("Y (m)")
ax.set_zlabel("Z (m)")
# ax.set_title("Pelota de golf rugosa vs lisa")
plt.savefig("golf_vssin.png", dpi=300)
plt.show()

# Vista desde un punto de vista superior del plano xy
plt.figure()
plt.plot(tray_hook[:,0], tray_hook[:,1], label='hook -- rugosa')
plt.plot(tray_slice[:,0], tray_slice[:,1], label='slice -- rugosa')
plt.plot(tray_no[:,0], tray_no[:,1], label='no Magnus -- rugosa')
plt.plot(tray_hookl[:,0], tray_hookl[:,1], label='hook -- lisa')
plt.plot(tray_slicel[:,0], tray_slicel[:,1], label='slice -- lisa')
plt.plot(tray_nol[:,0], tray_nol[:,1], label='no Magnus -- lisa')
plt.xlabel("X (m)")
plt.ylabel("Y (m)")
# plt.legend()
plt.savefig("golf_lateralXYsin.png", dpi=300, bbox_inches="tight")
# plt.title("Vista superior (desviación lateral)")
plt.show()

# Vista desde un punto de vista superior del plano xz
plt.figure()
plt.plot(tray_hook[:,0], tray_hook[:,2], label='hook -- rugosa')
plt.plot(tray_slice[:,0], tray_slice[:,2], label='slice -- rugosa')
plt.plot(tray_no[:,0], tray_no[:,2], label='no Magnus -- rugosa')
plt.plot(tray_hookl[:,0], tray_hookl[:,2], label='hook -- lisa')
plt.plot(tray_slicel[:,0], tray_slicel[:,2], label='slice -- lisa')
plt.plot(tray_nol[:,0], tray_nol[:,2], label='no Magnus -- lisa')
plt.xlabel("X (m)")
plt.ylabel("Z (m)")
# plt.legend()
# plt.title("Vista lateral (XZ)")
plt.savefig("golf_lateralXZsin.png", dpi=300, bbox_inches="tight")
plt.show()

# Vista desde un punto de vista superior del plano yz
plt.figure()
plt.plot(tray_hook[:,1], tray_hook[:,2], label='hook -- rugosa')
plt.plot(tray_slice[:,1], tray_slice[:,2], label='slice -- rugosa')
plt.plot(tray_no[:,1], tray_no[:,2], label='no Magnus -- rugosa')
plt.plot(tray_hookl[:,1], tray_hookl[:,2], label='hook -- lisa')
plt.plot(tray_slicel[:,1], tray_slicel[:,2], label='slice -- lisa')
plt.plot(tray_nol[:,1], tray_nol[:,2], label='no Magnus -- lisa')
plt.xlabel("Y (m)")
plt.ylabel("Z (m)")
# plt.legend()
# plt.title("Vista lateral (YZ)")
plt.savefig("golf_lateralYZsin.png", dpi=300, bbox_inches="tight")
plt.show()