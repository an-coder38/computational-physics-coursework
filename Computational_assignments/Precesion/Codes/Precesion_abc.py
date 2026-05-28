# Se parte del código subido a moodle
from math import pi, sqrt
from numpy import zeros
from vpython import sphere, rate, color, canvas, vector
import matplotlib.pyplot as plt


# Definición de parámetros de la simulación
GMSol = 4 * pi**2        # G * Masa del Sol en Unidades Astronómicas
dt = 0.001              # Paso temporal
T = 76                   # Periodo en AU
tiempoTotal = 1 * T     # Tiempo total del cálculo
NN = tiempoTotal / dt   # Número de pasos
N = int(NN)


# Creación de vectores de posición y velocidad
x =  zeros(N, float)
y =  zeros(N, float)
r =  zeros(N, float)    # Radio vector del planeta
vx = zeros(N, float)
vy = zeros(N, float)
t  = zeros(N, float)
v2  = zeros(N, float) #velocidad al cuadrado
U = zeros(N, float)
K = zeros(N, float)
E = zeros(N, float)
L = zeros(N, float)

# Condiciones Iniciales expresadas en unidades astronómicas y años.
t[0] = 0
x[0] = 0.59
y[0] = 0
r[0] = sqrt(x[0]**2 + y[0]**2)
vx[0] = 0
vy[0] = 11.47250862 # con este valor sale: 75.99899999999883 0.5899997651358467 1.8773405232166973e-05
#prueba y error hasta que solo de una vuelta y que además su última posición vuelva a ser la original

# energías iniciales y momento angular inicial todo por unidad de masa (supusimos m=0)
v2[0] = vx[0]**2+vy[0]**2
K[0] = 0.5*v2[0]
U[0] = -GMSol/r[0]
E[0] = K[0]+U[0]
L[0] = x[0]*vy[0]-y[0]*vx[0]

# Bucle principal utilizando método de Euler-Cromer
for i in range(N-1):
    vx[i+1] = vx[i] - GMSol * (x[i]/ r[i]**3) * dt
    vy[i+1] = vy[i] - GMSol * (y[i]/ r[i]**3) * dt
    x[i+1]  = x[i] + vx[i+1] * dt
    y[i+1]  = y[i] + vy[i+1] * dt
    r[i+1] = sqrt(x[i+1]**2 + y[i+1]**2)
    t[i+1] = t[i] + dt
    v2[i+1] = vx[i+1]**2+vy[i+1]**2

# obtener las energías y el momento angular
    K[i+1] = 0.5*v2[i+1]
    U[i+1] = -GMSol/r[i+1]
    E[i+1] = K[i+1]+U[i+1]
    L[i+1] = x[i+1]*vy[i+1]-y[i+1]*vx[i+1]

    if y[i+1]<0.1 and y[i+1]>=0 and x[i+1]>0:
        print(t[i+1], x[i+1], y[i+1])

print("Afelio:",max(r),"UA")
print("Velocidad máxima:",sqrt(max(v2)),"UA/año")


# Gráficas

plt.figure()
plt.plot(t,K,label="Ec/m")
plt.plot(t,U,label="Ep/m")
plt.plot(t,E,label="Et/m")
plt.legend()
plt.xlabel("Tiempo (años)")
plt.ylabel("E/m")
plt.savefig("energia_halley.png", dpi=300)
plt.show()

plt.figure()
plt.plot(t,L)
plt.xlabel("Tiempo (años)")
plt.ylabel("L/m")
plt.savefig("momento_halley.png", dpi=300)
plt.show()


# b) Plutón
print ('---------Apartado b----------')
# en vez de crear una función para ambos creo dos bucles diferentes y en este no coloco lo de
# mostrar los valores en la terminal ya que uso la velocidad exacta de la fórmula de la velocidad

#escribo como tiempo de simulación el del periodo de Plutón para comparar la forma de las energías con las de Halley
TP = 248.26 # con un tiempo superior al periodo de Plutón (en el dibujo de la órbita sale la órbita entera), puedo obtener el valor del periodo a partir de a
tiempoTotalP = 1 * TP     # Tiempo total del cálculo
NNP = tiempoTotalP / dt   # Número de pasos
NP = int(NNP)
# Creación de vectores de posición y velocidad para Plutón
xP =  zeros(NP, float)
yP =  zeros(NP, float)
rP =  zeros(NP, float)    # Radio vector del planeta
vxP = zeros(NP, float)
vyP = zeros(NP, float)
tP  = zeros(NP, float)
v2P  = zeros(NP, float) #velocidad al cuadrado
UP = zeros(NP, float)
KP = zeros(NP, float)
EP = zeros(NP, float)
LP = zeros(NP, float)


rpP = 29.7
raP = 49.3
aP = (rpP+raP)/2
vpP = sqrt(GMSol*(2/rpP - 1/aP)) #obtener la velocidad con la fórmula


tP[0] = 0
xP[0] = rpP
yP[0] = 0
rP[0] = sqrt(xP[0]**2 + yP[0]**2)
vxP[0] = 0
vyP[0] = vpP


# energías iniciales y momento angular inicial todo por unidad de masa (supusimos m=0)
v2P[0] = vxP[0]**2+vyP[0]**2
KP[0] = 0.5*v2P[0]
UP[0] = -GMSol/rP[0]
EP[0] = KP[0]+UP[0]
LP[0] = xP[0]*vyP[0]-yP[0]*vxP[0]

# Bucle principal utilizando método de Euler-Cromer
for i in range(NP-1):
    vxP[i+1] = vxP[i] - GMSol * (xP[i]/ rP[i]**3) * dt
    vyP[i+1] = vyP[i] - GMSol * (yP[i]/ rP[i]**3) * dt
    xP[i+1]  = xP[i] + vxP[i+1] * dt
    yP[i+1]  = yP[i] + vyP[i+1] * dt
    rP[i+1] = sqrt(xP[i+1]**2 + yP[i+1]**2)
    tP[i+1] = tP[i] + dt
    v2P[i+1] = vxP[i+1]**2+vyP[i+1]**2

# obtener las energías y el momento angular
    KP[i+1] = 0.5*v2P[i+1]
    UP[i+1] = -GMSol/rP[i+1]
    EP[i+1] = KP[i+1]+UP[i+1]
    LP[i+1] = xP[i+1]*vyP[i+1]-yP[i+1]*vxP[i+1]


# Resultados: r--Haley; rP--Plutón

print("Halley afelio:",max(r),"UA")
print("Pluton afelio:",max(rP),"UA")
print("Halley perihelio:",min(r),"UA")
print("Pluton perihelio:",min(rP),"UA")
print("Halley a:", 0.5*(min(r)+max(r)),"UA")
print("Pluton a:",0.5*(min(rP)+max(rP)),"UA")
print("Halley excentricidad:", (max(r)-min(r))/(min(r)+max(r)),"UA")
print("Pluton excentricidad:",(max(rP)-min(rP))/(min(rP)+max(rP)),"UA")
print("Halley T:", (0.5*(min(r)+max(r)))**1.5,"años") # 3LK, T^2=a^3
print("Pluton T:", (0.5*(min(rP)+max(rP)))**1.5,"años") # 3LK, T^2=a^3
#en nuestras unidades, mu=4*pi^2
print("Halley Et:", -0.5*4*(pi**2)/(0.5*(min(r)+max(r))),"años") # Et=-mu/2a
print("Pluton Et:", -0.5*4*(pi**2)/(0.5*(min(rP)+max(rP))),"años") # Et=-mu/2a

#Comparación visual de las órbitas
plt.figure(figsize=(7,7))
plt.plot(x,y,label="Cometa Halley")
plt.plot(xP,yP,label="Plutón")
plt.scatter(0,0,label="Sol")
plt.xlabel("x (UA)")
plt.ylabel("y (UA)")
plt.axis("equal")
plt.legend()
plt.savefig("halley_vs_pluton.png", dpi=300)
plt.show()

# Gráficas de plutón

plt.figure()
plt.plot(tP,KP,label="Ec/m")
plt.plot(tP,UP,label="Ep/m")
plt.plot(tP,EP,label="Et/m")
plt.legend()
plt.xlabel("Tiempo (años)")
plt.ylabel("E/m")
plt.savefig("energia_pluton.png", dpi=300)
plt.show()

plt.figure()
plt.plot(tP,LP)
plt.xlabel("Tiempo (años)")
plt.ylabel("L/m")
plt.savefig("momento_pluton.png", dpi=300)
plt.show()

# #Representación del movimiento planetario.
# canvas(title = "Órbita de la Tierra", x = 500, y = 100, width = 700, height = 700, center = vector(0, 0, 0.4), 
#         forward = vector(0, 0, -1), background = vector(0, 0, 0))
        
# sol = sphere(pos = vector(0.0, 0, 0), color = color.yellow, radius = 0.2)
# planeta = sphere(pos = vector(x[0], y[0], 0), color = color.cyan, radius = 0.05)

# #Bucle que representa el movimiento y nos proporciona el tiempo transcurrido.
# for i in range(N-1):
#     rate(50)
#     px = x[i]
#     py = y[i]
#     sphere(pos = vector(px, py, 0), color = color.red, radius = 0.01)
#     planeta.pos = vector(px, py, 0)
#     print(t[i])
    
