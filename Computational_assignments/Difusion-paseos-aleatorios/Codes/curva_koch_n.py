import matplotlib.pyplot as plt
import numpy as np

def koch_segment(p1, p2, n):
    # obtiene los segmentos de las curvas de Koch
    if n == 0:
        return [p1, p2]
    
    p1 = np.array(p1)
    p2 = np.array(p2)
    
    # Dividir el segmento en tres partes
    pA = p1 + (p2 - p1) / 3
    pB = p1 + 2 * (p2 - p1) / 3
    
    # Calcular el punto del pico (triángulo equilátero)
    angle = np.pi / 3  # 60 grados
    direction = pB - pA
    rotation_matrix = np.array([
        [np.cos(angle), -np.sin(angle)],
        [np.sin(angle),  np.cos(angle)]
    ])
    pC = pA + rotation_matrix.dot(direction)
    
    # Recursión en los 4 segmentos
    seg1 = koch_segment(p1, pA, n-1)
    seg2 = koch_segment(pA, pC, n-1)
    seg3 = koch_segment(pC, pB, n-1)
    seg4 = koch_segment(pB, p2, n-1)
    
    # Unir segmentos evitando puntos duplicados
    return seg1[:-1] + seg2[:-1] + seg3[:-1] + seg4

def dibujar_koch(n):
    # dibujar la curva de koch de orden n, para obtener la curva se utiliza la función koch_segment ya definida
    p1 = [0, 0]
    p2 = [1, 0]
    
    puntos = koch_segment(p1, p2, n)
    
    x, y = zip(*puntos)
    
    plt.figure(figsize=(8, 3))
    plt.plot(x, y)
    plt.title(f"Curva de Koch (orden {n})")
    plt.axis('equal')
    plt.axis('off')
    plt.show()
    return 0

def longitud_curva(puntos):
    # calcular la longitud de la curva de dimensión n
    puntos = np.array(puntos)
    diffs = np.diff(puntos, axis=0)
    distancias = np.linalg.norm(diffs, axis=1)
    return np.sum(distancias)

def area_poligono(n):
    # se van sumando las áreas de los polígonos que se crean
    A0 = np.sqrt(3)/4  # área triángulo inicial
    A = A0
    for k in range(1, n+1):
        N = 3 * 4**(k-1)
        a_tri = A0 / 9**k
        A += N * a_tri
    return A

def area_teorica_limite():
    A0 = np.sqrt(3)/4  # área del triángulo inicial (lado 1)
    return (8/5) * A0

def box_counting(puntos, epsilons):
    #método box_counting para posteriormente poder obtener la dimensión fractal de las curvas de koch
    puntos = np.array(puntos)
    
    counts = []
    
    for eps in epsilons:
        # discretizar espacio
        grid = np.floor(puntos / eps)
        # contar celdas únicas
        unique_boxes = np.unique(grid, axis=0)
        counts.append(len(unique_boxes))
    
    return counts

def dimension_fractal(puntos):
    epsilons = np.logspace(-3, -1, 10)
    counts = box_counting(puntos, epsilons)
    
    log_eps = np.log(1/epsilons)
    log_counts = np.log(counts)
    
    # ajuste lineal
    coef = np.polyfit(log_eps, log_counts, 1)
    return coef[0]

def estudio_koch(n):
    #función que implementa el resto de funciones ya definidas para estudiar la curva de koch de orden n
    
    dibujar_koch(n)

    puntos = koch_segment([0,0], [1,0], n)
    
    L = longitud_curva(puntos)
    longitud_teorica = (4/3)**n
    error_L = abs(L - longitud_teorica)
    D = dimension_fractal(puntos)
    error_D = abs(D - np.log(4)/np.log(3)) #log(4)/log(3) es la D real
    
    A = area_poligono(n)
    
    print(f"--------------Orden {n}--------------")
    print(f'Curva dibujada de orden: {n}')
    print(f"Longitud: {L}")
    print(f"Error en la longitud: {error_L}")    
    print(f"Área (copo): {A}")
    print(f"Dimensión fractal (numérica): {D}")
    print(f"Error en la dimensión fractal: {error_D}")

def grafica_longitud(n_max):
    ns = []
    Ls = []
    
    for n in range(n_max+1):
        puntos = koch_segment([0,0], [1,0], n)
        Ls.append(longitud_curva(puntos))
        ns.append(n)
    
    plt.plot(ns, Ls, marker='o')
    plt.title("Crecimiento de la longitud")
    plt.xlabel("n")
    plt.ylabel("Longitud")    
    plt.savefig("long_koch.png", dpi=300)
    plt.show()

def grafica_area(n_max):
    ns = []
    As = []
    
    for n in range(n_max+1):
        As.append(area_poligono(n))
        ns.append(n)
    
    plt.plot(ns, As, marker='o')
    plt.axhline(area_teorica_limite(), linestyle='--')
    plt.title("Convergencia del área")
    plt.xlabel("n")
    plt.ylabel("Área")
    plt.savefig("area_koch.png", dpi=300)
    plt.show()

def grafica_box_counting(puntos):
    epsilons = np.logspace(-3, -1, 10)
    counts = box_counting(puntos, epsilons)
    
    log_eps = np.log(1/epsilons)
    log_counts = np.log(counts)
    
    plt.scatter(log_eps, log_counts)
    
    coef = np.polyfit(log_eps, log_counts, 1)
    plt.plot(log_eps, np.polyval(coef, log_eps))
    
    plt.xlabel("log(1/ε)")
    plt.ylabel("log(N(ε))")
    plt.title("Estimación de dimensión fractal")
    plt.show()

def grafica_dimension(n_max):
    ns = []
    Ds = []
    
    D_teorica = np.log(4)/np.log(3)
    
    for n in range(n_max+1):
        puntos = koch_segment([0,0], [1,0], n)
        D = dimension_fractal(puntos)
        
        ns.append(n)
        Ds.append(D)
    
    plt.plot(ns, Ds, marker='o', label="Dimensión numérica")
    plt.axhline(D_teorica, linestyle='--', label="Dimensión teórica")
    plt.xlabel("n")
    plt.ylabel("Dimensión fractal")
    plt.title("Convergencia de la dimensión fractal")
    plt.savefig("dim_fractal.png", dpi=300)
    plt.legend()
    plt.show()



n_max = 10 # a partir de n=10, el tiempo de ejecución aumenta drásticamente
for n in range(n_max+1):    
    estudio_koch(n)
    A = area_poligono(n)
    print(n, A)

grafica_longitud(n_max)
grafica_area(n_max)
grafica_dimension(n_max)
grafica_area(20)

