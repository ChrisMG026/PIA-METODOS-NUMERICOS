# Manual técnico — PIA de Métodos Numéricos

## Propósito y ejecución

El programa reconstruye los parámetros de los documentos A–G del expediente del PIA a partir de cuatro matrículas de siete dígitos. El punto de entrada es `PIA.py`; se ejecuta con `python PIA.py`. Requiere Python 3.10 o posterior con Tkinter. No requiere dependencias externas.

La aplicación ordena las matrículas de menor a mayor antes de asignarles M1–M4. Para cada matrícula calcula:

- \(M_i\): valor entero de sus siete dígitos.
- \(N_i=M_i/10^7\).
- \(P_i=(1/7)\sum_{k=1}^7d_{ik}\).
- \(S_i=\sum_{k=1}^7k\,d_{ik}\).

## Métodos implementados

Todos los métodos numéricos empleados se implementan en `PIA.py`:

- Series de Mercator transformada para \(\ln(x)\) y de Taylor para \(\exp(x)\), con tolerancia de término \(10^{-12}\).
- Interpolación de Lagrange.
- Newton–Raphson para C y \(A_r\), con historial de aproximaciones.
- Mínimos cuadrados para las pendientes globales \(C_r\), \(k\) y \(\chi\).
- Eliminación de Gauss–Jordan con pivoteo parcial para los sistemas de C.
- Regla del trapecio compuesta, Simpson 1/3 y Simpson 3/8.
- Bisección para \(n\), con tolerancia \(10^{-10}\), intervalo abierto \(1<n<3\) y registro de iteraciones.
- Euler hacia adelante para las EDO de D, E y F. En E y F se conserva el historial completo \((t,y)\).
- Aislamiento por intervalos y bisección de las raíces reales de la cúbica de saturación.

## Reconstrucción de los documentos

### A — Estado operativo

Se calculan las variables \(a,\ldots,i,A,B,C\) con las expresiones y tablas del expediente. Para \(A\) se integra \(2t^2+S_2t+1\) entre 0 y \(P_1\); para \(B\) se interpola la tabla en \(N_2\); C usa Newton–Raphson para la cúbica recuperada. Se muestran las iteraciones de C.

### B — Respuesta bajo carga

Se calculan \(A_r\) como la menor raíz positiva de la cúbica del expediente, \(B_r\) por la integral indicada, \(C_r\) como pendiente de mínimos cuadrados sobre la tabla recuperada y \(D_r\) por Lagrange en \(x=2.5\).

### C — Saturación

Se calculan \(a_s=R(2)\), \(b_s=x+y+z\) mediante el sistema 3×3 recuperado, \(c_s\) con trapecio compuesto en \([0,4]\), y \(d_s=10T_m+3T_c\). \(T_m\) se interpola en 410 y \(T_c\) se evalúa en el parámetro \(a\) del Documento A.

### D — Balance acumulado

La pendiente global \(k\) usa los puntos \((1,2.0),(3,5.3),(5,8.9),(7,12.2)\). También se calculan \(m\) con Simpson 1/3, \(n\) con bisección en \(1<n<3\), \(p\) con Euler desde \(y(0)=1\) y paso 0.1 hasta \(t=1.5\), y \(q\) con Lagrange en \(r=4.5\).

### E — Dinámica del rendimiento

Se elige \(r_\sigma=1.2\), la menor raíz positiva de \(r^2-5.2r+4.8=0\), y se calcula \(\sigma=r_\sigma/10\). La ecuación se integra de 0 a 3 con Euler y paso 0.1; \(Y(0)=R(1.8)+x\).

**Supuesto declarado:** el expediente no da un valor independiente para \(\omega\). Se conserva la interpretación del código de trabajo, \(\omega=\kappa=L(3)\), así que \(\kappa/\omega=1\). Debe confirmarse este supuesto con el docente.

### F — Comportamiento térmico

Se calcula \(\chi\) por mínimos cuadrados, \(\psi=8\) como raíz mayor de \(r^2-11r+24=0\), y \(\xi=a+d+g\) del Documento A. Se propone \(\theta(0)=25\,^\circ\mathrm C\), como temperatura ambiente o de reposo antes de iniciar la carga de trabajo. Se integra hasta \(t=3\) con Euler y paso 0.1.

**Corrección explícita respecto del expediente:** la ecuación transcrita, \(\theta'=\chi\theta(\psi-\theta)-Y/\xi\), da divergencia para entradas típicas: \(Y/\xi\) supera la capacidad logística máxima \(\chi\psi^2/4\). La integración directa produce temperaturas no finitas antes de completar el intervalo, por lo que tampoco existe un valor numérico válido de \(\theta(1)\) para G.

El signo negativo de \(Y/\xi\) hace que ese término reduzca la temperatura en la ecuación transcrita; por ello, no se debe describir como un término que eleva la temperatura. La condición inicial de 25 °C representa el estado previo a la carga, no una conclusión sobre el signo del cambio posterior.

Para conservar la condición inicial propuesta y obtener una trayectoria finita, el programa normaliza la fuente:

\[
\theta'=\chi\theta(\psi-\theta)-\frac{Y}{\xi\lambda},\qquad
\lambda=\max\left(1,\frac{\max_t(Y(t)/\xi)}{0.8\,\chi\psi^2/4}\right).
\]

El factor \(\lambda\) limita la fuente al 80 % de la capacidad logística máxima, dejando un margen numérico del 20 %. Este factor y el margen son una decisión de modelado añadida por el equipo, no datos del PDF. Deben reportarse como desviación y validarse con el docente antes de entregar resultados como si provinieran literalmente del expediente.

### G — Índice de estabilidad

Se resuelve \(a_sU^3-b_sU^2+c_sU-d_s=0\). Si hay varias raíces reales positivas, \(U^*\) es la menor; se declara este criterio porque el expediente no especifica cuál solución operativa elegir. La integral de \(Y(t)\) en \([0,2]\) se calcula con Simpson 1/3 usando la trayectoria de E; \(\theta(1)\) se obtiene interpolando la trayectoria de F. Finalmente:

\[
S=\tfrac12\int_0^2Y(t)\,dt+\tfrac13\theta(1)+\ln(1+U^*).
\]

## Validación y sensibilidad

La pantalla **Análisis y gráficas** y el informe generado comparan Euler con \(h=0.1,0.05,0.025\), muestran las diferencias consecutivas de \(S\) y repiten Newton–Raphson para C con tolerancias \(10^{-6},10^{-8},10^{-10}\). Las gráficas se construyen directamente desde las trayectorias calculadas, no desde datos externos.

El menú **Operaciones detalladas** presenta sumatorias, iteraciones de Newton y bisección, matrices intermedias de Gauss–Jordan, resultados y supuestos. El botón **Exportar informe** guarda un HTML reproducible con tablas y gráficas SVG; se puede imprimir desde el navegador como PDF.

El HTML exportado incluye las matrículas de la instancia, resultados A–G, sensibilidad, validaciones numéricas, interpretación, supuestos y el contenido completo del archivo ejecutable. Imprimir ese HTML a PDF genera el informe técnico asociado a la instancia capturada.

## Limitaciones que deben reportarse

1. El expediente no especifica \(\omega\); el programa usa \(\omega=\kappa\).
2. La ecuación F requiere la normalización \(\lambda\) descrita para evitar divergencia. Es un supuesto añadido.
3. La condición inicial \(\theta(0)=25\,^\circ\mathrm C\) es propuesta por el equipo, no recuperada del archivo.
4. La selección de la menor raíz positiva para \(U^*\) es un criterio documentado, no una selección conservada en el expediente.
5. Las soluciones de Euler dependen de \(h\); el informe muestra una comparación, no una cota rigurosa del error global.
