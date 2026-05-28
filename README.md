# Computational Physics & Advanced Computing | UAM

Este repositorio contiene las simulaciones computacionales, métodos numéricos e informes científicos desarrollados para la asignatura de **Computación Avanzada** en el Grado en Física de la Universidad Autónoma de Madrid (UAM).

## Descripción General
El núcleo del repositorio consiste en el diseño e implementación de soluciones numéricas, algoritmos deterministas y métodos estocásticos aplicados a sistemas físicos complejos donde las aproximaciones analíticas son insuficientes. Los proyectos cubren un amplio espectro de la física fundamental: la transición al caos determinista, hasta la propagación de paquetes de ondas cuánticas.

## Estructura del Repositorio

El proyecto está estructurado para garantizar la reproducibilidad y separar el entorno de desarrollo del material académico definitivo:

* **`exercises/`**: Scripts breves y preámbulos algorítmicos resueltos durante las sesiones de clase. Cada carpeta numérica incluye sus respectivas hojas de instrucciones (`Guidelines/`) y los códigos de resolución (`Solution/`).
* **`Computational_assignments/`**: Proyectos de investigación principales. Cada módulo temático está aislado en su propia carpeta con la siguiente estructura interna:
  * `Guidelines/`: Enunciados oficiales, requerimientos y fundamentos teóricos de la práctica.
  * `Codes/`: Scripts en Python (`.py`) optimizados para la ejecución de las simulaciones, análisis de datos y ampliaciones avanzadas.
  * `Reports/`: Memorias científicas formales que incluyen el informe final compilado en **PDF**, el documento fuente en **LaTeX** (`.tex`) y el repositorio de imágenes y gráficas generadas.

---

## Módulos Temáticos Desarrollados

Los proyectos principales integrados en `Computational_assignments/` se organizan cronológicamente según su complejidad física y algorítmica:

### Cinemática y Aerodinámica de Proyectiles Reales (`Proyectiles`)
* **Dinámica en 2D (Béisbol):** Implementación del método de Euler para resolver ecuaciones diferenciales no lineales ordinarias. Integración de la fuerza de arrastre (drag) dependiente de la velocidad y el efecto de corrientes de viento lateral.
* **Dinámica en 3D (Golf):** Incorporación del efecto Magnus, analizando cómo el espín (rotación de la bola) y la rugosidad superficial alteran la sustentación, el alcance máximo, el ángulo óptimo de lanzamiento y la desviación lateral.

### Dinámica No Lineal y Caos Determinista (`Caos`)
* **Péndulo Forzado y Amortiguado:** Estudio de la transición al caos en sistemas dinámicos no lineales sometidos a una fuerza impulsora periódica externa.
* **Espacio de Fases y Bifurcaciones:** Construcción de secciones de Poincaré y diagramas de bifurcación mediante barridos paramétricos de la fuerza externa, analizando la ruta de duplicación de período (*period-doubling*).
* **Exponente de Lyapunov:** Estimación numérica cuantitativa del exponente de Lyapunov ($\lambda$) para caracterizar la sensibilidad extrema a las condiciones iniciales y evaluar cómo el coeficiente de amortiguamiento mitiga el régimen caótico.

### Mecánica Celeste y Perturbaciones Gravitatorias (`Precesion`)
* **Simulación Orbital:** Modelización de la dinámica del cometa Halley bajo el potencial gravitatorio del Sol y comparación geométrica/energética de su órbita altamente elíptica con la de Plutón.
* **Perturbaciones de 3 Cuerpos:** Evaluación del efecto gravitatorio de Júpiter sobre el cometa mediante un método de extrapolación lineal (escalado artificial de la masa joviana) para cuantificar la precesión inducida del perihelio.
* **Extensión Relativista:** Verificación cualitativa de las anomalías orbitales aplicadas a la precesión de Mercurio en la mecánica clásica frente a correcciones perturbativas.

### Procesos Estocásticos, Difusión y Geometría Fractal (`Difusion-paseos-aleatorios`)
* **Equivalencia Micro/Macro:** Demostración de la correspondencia estadística entre el paseo aleatorio discreto en dos dimensiones (crecimiento lineal del desplazamiento cuadrático medio) y la solución continua de la ecuación de difusión.
* **Evolución Entrópica:** Modelización de la entropía de un sistema de caminantes independientes bajo condiciones periódicas de frontera, cuantificando la convergencia irreversible hacia el estado de equilibrio homogéneo.
* **Sistemas Fractales:** Construcción geométrica de la curva de Koch y estimación numérica de su dimensión fractal no entera mediante el algoritmo de conteo por cajas (*box-counting*).

### Fenómenos Colectivos y Mecánica Estadística (`Ising`)
* **Transiciones de Fase:** Simulación de un ferromagneto bidimensional a través del Modelo de Ising 2D en una red cuadrada de espines utilizando el algoritmo de muestreo de Metropolis.
* **Histéresis Magnética:** Análisis del comportamiento de la magnetización por espín $M(H)$ ante ciclos de campos magnéticos externos variables, evaluando la influencia crítica del tiempo de relajación de Monte Carlo.
* **Termodinámica Crítica:** Caracterización de las discontinuidades en la energía interna $E(T)$ y la susceptibilidad magnética en las inmediaciones de la temperatura crítica de Onsager ($T_c \approx 2.26$).

### Problemas de Contorno y Métodos de Relajación (`Relajacion`)
* **Algoritmos de Relajación:** Implementación del método iterativo de Gauss-Seidel para la resolución de ecuaciones diferenciales sujetas a condiciones de contorno espaciales o temporales puras.
* **Estudio de Casos:**
  * *Dinámica Unidimensional:* Trayectoria de una partícula bajo aceleración constante con restricciones fijadas exclusivamente en los extremos temporales del movimiento.
  * *Electrostática:* Resolución de la ecuación de Laplace en coordenadas cilíndricas con simetría axial para determinar el perfil de potencial electrostático en un condensador cilíndrico tras su discretización espacial.

### Dinámica Cuántica y Difracción de Materia (`Cuantica`)
* **Ecuación de Schrödinger:** Resolución numérica de la ecuación de Schrödinger dependiente del tiempo para analizar la evolución espacial y temporal de un paquete de ondas gausiano.
* **Efecto Túnel:** Cuantificación continua de los coeficientes de transmisión y reflexión de la densidad de probabilidad al interactuar con barreras de potencial rectangulares e impenetrables.
* **Difracción Macromolecular:** Extensión bidimensional para simular el experimento de la doble rendija aplicado a la frontera de la física cuántica, modelando por métodos estocásticos (Monte Carlo) la difracción coherente de moléculas masivas de nanografeno (Hexabenzocoroneno, $C_{42}H_{18}$).

---

## Stack Técnico
* **Lenguaje:** Python
* **Gestor de Entornos:** Miniconda3 / Archivos de entorno específicos
* **Librerías Científicas:** `NumPy`, `SciPy`, `Matplotlib`
* **Entornos de Desarrollo:** VS Code / Spyder
* **Redacción Científica:** LaTeX (Overleaf)

## Autor
**Andrés López Serna** – Estudiante de Grado en Física, Universidad Autónoma de Madrid (UAM).