import numpy as np 
import matplotlib.pyplot as plt

pi=np.pi
sqrt = np.sqrt
atan2 = np.atan2

def filtrar_outliers(x, y, n_sigma=2.5):
    x = np.array(x)
    y = np.array(y)
    media = np.mean(y)
    std = np.std(y)
    mask = np.abs(y - media) < n_sigma*std
    
    return x[mask], y[mask]

def filtrar_outliers_mad(x, y, thresh=3.5):
    x = np.array(x)
    y = np.array(y)
    med = np.median(y)
    mad = np.median(np.abs(y - med))
    if mad == 0:
        return x, y

    z = 0.6745 * (y - med) / mad
    mask = np.abs(z) < thresh

    return x[mask], y[mask]

# funcion aceleracion Halley
def accel_halley(xh,yh,xj,yj): #MJ
    r = sqrt(xh**2 + yh**2)
    dx = xh-xj
    dy = yh-yj
    d = sqrt(dx**2 + dy**2) #distancia de halley a júpiter

    ax = -GMS*xh/r**3 - GMS*(MJ/MS)*dx/d**3
    ay = -GMS*yh/r**3 - GMS*(MJ/MS)*dy/d**3

    return ax,ay

# funcion aceleracion Jupiter
def accel_jupiter(xh,yh,xj,yj):
    r = sqrt(xj**2 + yj**2)
    dx = xj-xh
    dy = yj-yh
    d = sqrt(dx**2 + dy**2) #distancia de halley a júpiter

    ax = -GMS*xj/r**3 - GMS*(MH/MS)*dx/d**3
    ay = -GMS*yj/r**3 - GMS*(MH/MS)*dy/d**3

    return ax,ay


def sistema_precesion_ec (dt, T, t0, x0,y0,vx0,vy0,xj0,yj0,vxj0,vyj0): #m
    
    # MJ = 0.0009543*m
        
    N = int(T/dt)

    # arrays Halley
    x = np.zeros(N)
    y = np.zeros(N)
    vx = np.zeros(N)
    vy = np.zeros(N)
    t  = np.zeros(N, float)
    # r =  np.zeros(N, float)    # Radio vector del planeta
    # v2  = np.zeros(N, float) #velocidad al cuadrado
    # U = np.zeros(N, float)
    # K = np.zeros(N, float)
    # E = np.zeros(N, float)
    # L = np.zeros(N, float)

    # arrays Jupiter
    xj = np.zeros(N)
    yj = np.zeros(N)
    vxj = np.zeros(N)
    vyj = np.zeros(N)
    
    # condiciones iniciales Halley
    t[0]=t0
    x[0] = x0
    y[0] = y0
    vx[0] = vx0
    vy[0] = vy0
    
    # r[0] = sqrt(x[0]**2 + y[0]**2)
    # v2[0] = vx[0]**2+vy[0]**2
    # K[0] = 0.5*v2[0]
    # U[0] = -GMS/r[0]
    # E[0] = K[0]+U[0]
    # L[0] = x[0]*vy[0]-y[0]*vx[0]
    
    # condiciones iniciales Jupiter
    xj[0] = xj0
    yj[0] = yj0
    vxj[0] = vxj0
    vyj[0] = vyj0 #suponemos excentricidad de la órbita de júpiter igual a 0

    # listas perihelio
    perihelio_ang = []
    perihelio_t = []
    
    for i in range(N-1):

        # Halley (Euler-Cromer)
        ax_h,ay_h = accel_halley(x[i],y[i],xj[i],yj[i]) #MJ

        vx[i+1] = vx[i] + ax_h*dt
        vy[i+1] = vy[i] + ay_h*dt
        
        x[i+1] = x[i] + vx[i+1]*dt
        y[i+1] = y[i] + vy[i+1]*dt

        # Júpiter (Euler-Cromer)
        ax_j,ay_j = accel_jupiter(x[i],y[i],xj[i],yj[i])

        vxj[i+1] = vxj[i] + ax_j*dt
        vyj[i+1] = vyj[i] + ay_j*dt
        
        xj[i+1] = xj[i] + vxj[i+1]*dt
        yj[i+1] = yj[i] + vyj[i+1]*dt

    #     r[i+1] = sqrt(x[i+1]**2 + y[i+1]**2)
    #     t[i+1] = t[i] + dt
    #     v2[i+1] = vx[i+1]**2+vy[i+1]**2

    # # obtener las energías y el momento angular
    #     K[i+1] = 0.5*v2[i+1]
    #     U[i+1] = -GMS/r[i+1]
    #     E[i+1] = K[i+1]+U[i+1]
    #     L[i+1] = x[i+1]*vy[i+1]-y[i+1]*vx[i+1]


        # detectar perihelio del cometa Halley
        if i>2:

            r1 = sqrt(x[i-1]**2+y[i-1]**2)
            r2 = sqrt(x[i]**2+y[i]**2)
            r3 = sqrt(x[i+1]**2+y[i+1]**2)

            if r2<r1 and r2<r3:

                ang = atan2(y[i],x[i])

                perihelio_ang.append(ang)
                perihelio_t.append(i*dt)

    # # calcular precesion
    # if len(perihelio_ang)>1:

    #     precesion = perihelio_ang[-1]-perihelio_ang[0]

    #     print("Número de órbitas:",len(perihelio_ang))
    #     print("Precesión total:",precesion,"rad")
    #     print("Precesión por órbita:",precesion/(len(perihelio_ang)-1),"rad")

    
    return perihelio_ang,perihelio_t


# constantes
GMS = 4*pi**2
MS = 1
MH=0 # suponemos que el cometa no tiene masa

# condiciones iniciales
dt = 0.001
T = 500
t0=0
x0=0.59;y0=0;vx0=0;vy0=11.47250862 #CI Halley
aj = 5.2 #afelio júpiter


###########LO COMENTADO EN ESTA PARTE ES LO QUE HAY QUE CORRER PARA OBTENER EL .txt TARDA EN CORRER##########################
# m_values = np.linspace (1,10,80)
# fit = []

# for m in m_values:
    
#     MJ = 0.0009543*m
#     GMJ = GMS*MJ
#     xj0=5.2;yj0=0;vxj0=0;vyj0=sqrt(GMS/aj) #(1+GMJ/GMS)*GMS/aj) #CI  Júpiter
    
#     # t,x,y,vx,vy,r,v2,K,U,E,L,xj,yj,vxj,vyj,perihelio_ang,perihelio_t = sistema (dt, T, t0, x0,y0,vx0,vy0,xj0,yj0,vxj0,vyj0)
#     perihelio_ang,perihelio_t = sistema_precesion_ec(dt, T, t0, x0, y0, vx0, vy0, xj0, yj0, vxj0, vyj0)
#     pendiente, ordenada = np.polyfit(perihelio_t, perihelio_ang, 1)
#     fit.append(pendiente)
#     print(f'{m}')

# # juntar columnas
# datos = np.column_stack((m_values, fit))
# # guardar archivo
# np.savetxt("datos_m_fit_ec.txt", datos, header="m_values fit")


# cargar datos del archivo
data = np.loadtxt("datos_m_fit_ec.txt")
m_values = data[:,0]
fit = data[:,1]


# convertir a arrays
m_values = np.array(m_values)
fit = np.array(fit)
# filtrar inestabilidades
m_filtrado, fit_filtrado = filtrar_outliers_mad(m_values, fit)

print("EC--Datos originales:", len(m_values))
print("EC--Datos usados en el ajuste:", len(m_filtrado))

# ajuste final
pendiente_fin, ordenada_fin = np.polyfit(m_filtrado, fit_filtrado, 1)
# linea de ajuste
m_line = np.linspace(min(m_filtrado), max(m_filtrado), 10)
fit_line = pendiente_fin*m_line + ordenada_fin

# muestro todos los datos y los distingo de los filtrados
plt.figure()
plt.scatter(m_values, fit, alpha=0.25, label="Datos originales")
plt.scatter(m_filtrado, fit_filtrado, label="Datos usados")
plt.plot(m_line, fit_line, label="Ajuste lineal")
plt.xlabel("MJ usado / MJ real")
plt.ylabel("Precesión")
plt.title("EC--Extrapolación de la precesión del cometa Halley")
plt.legend()
plt.savefig("prec_halley_todo_ec.png", dpi=300)
plt.show()

# solo muestro los datos filtrados
plt.figure()
plt.scatter(m_filtrado, fit_filtrado, label="Datos usados")
plt.plot(m_line, fit_line, label="Ajuste lineal")
plt.xlabel("MJ usado / MJ real")
plt.ylabel("Precesión")
plt.title("EC--Extrapolación de la precesión del cometa Halley filtrados")
plt.legend()
plt.savefig("prec_halley_filtrado_ec.png", dpi=300)
plt.show()
print('EC--Precesión cometa Halley filtrado:', pendiente_fin*1 + ordenada_fin)


# quiero filtrar nuevamente únicamente para valores de m entre 1.3 y 6, donde parece que sí se mantiene la tendencia
# luego lo volveré a ajustar
mask_range = (m_filtrado > 1.3) & (m_filtrado < 6)

# aplicar máscara
m_filtrado = m_filtrado[mask_range]
fit_filtrado = fit_filtrado[mask_range]

# ajuste final
pendiente_fin, ordenada_fin = np.polyfit(m_filtrado, fit_filtrado, 1)
# linea de ajuste
m_line = np.linspace(min(m_filtrado), max(m_filtrado), 10)
fit_line = pendiente_fin*m_line + ordenada_fin

# solo muestro los datos filtrados
plt.figure()
plt.scatter(m_filtrado, fit_filtrado, label="Datos usados")
plt.plot(m_line, fit_line, label="Ajuste lineal")
plt.xlabel("MJ usado / MJ real")
plt.ylabel("Precesión")
plt.title("EC--Extrapolación de la precesión del cometa Halley: rango m = (1.3, 6)")
plt.legend()
plt.savefig("prec_halley_rango_ec.png", dpi=300)
plt.show()
print('EC--Precesión cometa Halley filtrado en el rango m = (1.3, 6):', pendiente_fin*1 + ordenada_fin,'rad/año')

precesion_rad_anio = pendiente_fin * 1.0 + ordenada_fin # extrapolar a MJ/MJ_real = 1
precesion_arcsec_siglo = precesion_rad_anio * (180/pi) * 3600 * 100
print('EC--Precesión cometa Halley filtrado en el rango m = (1.3, 6):', precesion_arcsec_siglo,'arcseg/siglo')
print('EC--Precesión cometa Halley filtrado en el rango m = (1.3, 6):', precesion_arcsec_siglo/3600,'grados/siglo')
print('pendiente(fit)', pendiente_fin)
print('ordenada(fit)', ordenada_fin)