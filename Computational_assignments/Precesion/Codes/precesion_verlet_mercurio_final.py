# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 09:36:05 2026

@author: jfcte
"""

import numpy as np
import matplotlib.pyplot as plt

pi = np.pi
sqrt = np.sqrt
atan2=np.atan2

# =========================
# CONSTANTES
# =========================

GMS = 4*pi**2      # AU^3 / yr^2
MS = 1
MM = 1.66e-7       # masa Mercurio
MJ_real = 0.0009543

# elementos orbitales Mercurio
a_mer = 0.387
e_mer = 0.2056
r_peri = a_mer*(1-e_mer)

# velocidad exacta en perihelio
v_peri = sqrt(GMS*(2/r_peri - 1/a_mer))

# Júpiter
a_jup = 5.2

def accel_mercury(xm, ym, xj, yj, MJ):

    r = sqrt(xm**2 + ym**2)

    dx = xm - xj
    dy = ym - yj
    d = sqrt(dx**2 + dy**2)

    ax = -GMS*xm/r**3 - GMS*(MJ/MS)*dx/d**3
    ay = -GMS*ym/r**3 - GMS*(MJ/MS)*dy/d**3

    return ax, ay


def accel_jupiter(xj, yj):

    r = sqrt(xj**2 + yj**2)

    ax = -GMS*xj/r**3
    ay = -GMS*yj/r**3

    return ax, ay


def sim_mercury(dt, T, MJ):

    N = int(T/dt)

    x = np.zeros(N)
    y = np.zeros(N)
    vx = np.zeros(N)
    vy = np.zeros(N)

    xj = np.zeros(N)
    yj = np.zeros(N)
    vxj = np.zeros(N)
    vyj = np.zeros(N)

    t = np.zeros(N)

    # condiciones iniciales Mercurio
    x[0] = r_peri
    y[0] = 0
    vx[0] = 0
    vy[0] = v_peri

    # condiciones iniciales Jupiter
    xj[0] = a_jup
    yj[0] = 0
    vxj[0] = 0
    vyj[0] = sqrt(GMS/a_jup)

    ax, ay = accel_mercury(x[0], y[0], xj[0], yj[0], MJ)
    axj, ayj = accel_jupiter(xj[0], yj[0])

    peri_ang = []
    peri_t = []

    for i in range(N-1):

        # posiciones
        x[i+1] = x[i] + vx[i]*dt + 0.5*ax*dt**2
        y[i+1] = y[i] + vy[i]*dt + 0.5*ay*dt**2

        xj[i+1] = xj[i] + vxj[i]*dt + 0.5*axj*dt**2
        yj[i+1] = yj[i] + vyj[i]*dt + 0.5*ayj*dt**2

        # nuevas aceleraciones
        ax_new, ay_new = accel_mercury(x[i+1], y[i+1], xj[i+1], yj[i+1], MJ)
        axj_new, ayj_new = accel_jupiter(xj[i+1], yj[i+1])

        # velocidades
        vx[i+1] = vx[i] + 0.5*(ax + ax_new)*dt
        vy[i+1] = vy[i] + 0.5*(ay + ay_new)*dt

        vxj[i+1] = vxj[i] + 0.5*(axj + axj_new)*dt
        vyj[i+1] = vyj[i] + 0.5*(ayj + ayj_new)*dt

        ax, ay = ax_new, ay_new
        axj, ayj = axj_new, ayj_new

        t[i+1] = t[i] + dt

        # detección de perihelio
        if i > 2:

            r1 = sqrt(x[i-1]**2 + y[i-1]**2)
            r2 = sqrt(x[i]**2 + y[i]**2)
            r3 = sqrt(x[i+1]**2 + y[i+1]**2)

            if r2 < r1 and r2 < r3:

                if len(peri_t) == 0 or (t[i]-peri_t[-1]) > 0.2:

                    peri_ang.append(atan2(y[i], x[i]))
                    peri_t.append(t[i])

    peri_ang = np.unwrap(np.array(peri_ang)) #unwrap elimina saltos bruscos de +-2pi
    peri_t = np.array(peri_t)

    return peri_ang, peri_t


# Extrapolación

m_values = np.linspace(2,10,20)
precesiones = []

dt = 1e-5
T = 100

###########LO COMENTADO EN ESTA PARTE ES LO QUE HAY QUE CORRER PARA OBTENER EL .txt TARDA EN CORRER##########################

# for m in m_values:

#     MJ = MJ_real*m

#     peri_ang, peri_t = sim_mercury(dt,T,MJ)

#     pendiente,_ = np.polyfit(peri_t, peri_ang,1)

#     precesiones.append(pendiente)

#     print("m =",m," pendiente =",pendiente)

# m_values = np.array(m_values)
# precesiones = np.array(precesiones)

# # Guardar datos de la extrapolación
# datos = np.column_stack((m_values, precesiones))
# np.savetxt("precesion_extrapolada_mercurio_end.txt", datos, header="MJ/MJ_real  pendiente[rad/yr]")

data = np.loadtxt("precesion_extrapolada_mercurio_end.txt")
m_values = data[:,0]
precesiones = data[:,1]

# ajuste lineal
pendiente_fit,ordenada = np.polyfit(m_values,precesiones,1)

precesion_rad_anio = pendiente_fit*1 + ordenada
precesion_arcsec_siglo = precesion_rad_anio*(180/pi)*3600*100

print("\nPrecesión extrapolada:",precesion_arcsec_siglo,"arcsec/siglo")


# Fraficar extrapolación para mercurio

plt.figure()

plt.scatter(m_values,precesiones,label="Simulaciones")

plt.plot(
    m_values,
    pendiente_fit*m_values+ordenada,
    'r--',
    label="Ajuste lineal"
)

plt.xlabel("MJ / MJ_real")
plt.ylabel("Precesión [rad/yr]")
plt.title("Extrapolación de la precesión de Mercurio")
plt.savefig("prec_mercurio_v.png", dpi=300)

plt.legend()
plt.grid()

plt.show()

# Para comprobar algunos resultados

# m_0 = 0
# MJ_0 = MJ_real*m_0

# peri_ang, peri_t = sim_mercury(dt,T,MJ_0)

# pendiente,_ = np.polyfit(peri_t, peri_ang,1)

# precesion_arcsec_siglo = pendiente*(180/pi)*3600*100

# print("m =",m_0," pendiente =",precesion_arcsec_siglo)

# m_1 = 1
# MJ_1 = MJ_real*m_1
# peri_ang, peri_t = sim_mercury(dt,T,MJ_1)

# pendiente,_ = np.polyfit(peri_t, peri_ang,1)

# precesion_arcsec_siglo = pendiente*(180/pi)*3600*100

# print("m =",m_1," pendiente =",precesion_arcsec_siglo)



