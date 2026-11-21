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
- Interpolación de Lagrange para nodos no equidistantes y Newton hacia adelante para tablas equidistantes.
- Newton–Raphson para C y aislamiento/bisección de raíces reales para los polinomios, con criterio explícito de selección de la menor raíz positiva para \(A_r\).
- Mínimos cuadrados para las pendientes globales \(C_r\), \(k\) y \(\chi\).
- Eliminación de Gauss–Jordan con pivoteo parcial para los sistemas de C.
- Simpson 1/3 compuesto (sobre funciones y sobre tablas) y Simpson 3/8.
- Bisección para \(n\), con tolerancia \(10^{-10}\), intervalo abierto \(1<n<3\) y registro de iteraciones.
- Runge-Kutta de cuarto orden para las EDO de D, E y F. En E y F se conserva el historial \((t,y)\); F conserva el historial parcial si la integración deja de ser finita.
- Aislamiento por intervalos y bisección de las raíces reales de la cúbica de saturación.

## Reconstrucción de los documentos

### A — Estado operativo

Se calculan las variables \(a,\ldots,i,A,B,C\) con las expresiones y tablas del expediente. Para \(A\) se integra \(2t^2+S_2t+1\) entre 0 y \(P_1\); para \(B\) se interpola la tabla en \(N_2\); C usa Newton–Raphson para la cúbica recuperada. Se muestran las iteraciones de C.

### B — Respuesta bajo carga

Se calculan \(A_r\) como la menor raíz positiva de la cúbica del expediente, \(B_r\) por la integral indicada, \(C_r\) como pendiente de mínimos cuadrados sobre la tabla recuperada y \(D_r\) por Newton hacia adelante en \(x=2.5\), pues la tabla tiene espaciamiento uniforme (el resultado coincide con Lagrange).

### C — Saturación

Se calculan \(a_s=R(2)\), \(b_s=x+y+z\) mediante el sistema 3×3 recuperado, \(c_s\) con Simpson 1/3 compuesto en \([0,4]\), y \(d_s=10T_m+3T_c\). \(T_m\) se interpola en 410 y \(T_c\) se evalúa en el parámetro \(a\) del Documento A.

Para \(c_s\) se usa la tabla \(\ell(t)=(3.0,4.4,6.9,10.1,15.0)\) en \(t=0,1,2,3,4\). Hay \(n=4\) intervalos (par, por lo que Simpson 1/3 compuesto es aplicable) y paso \(h=1\), con pesos \(1,4,2,4,1\):

\[
\sum_{\text{extremos}}=\ell_0+\ell_4=18,\qquad
\sum_{\text{impares}}=\ell_1+\ell_3=14.5,\qquad
\sum_{\text{pares}}=\ell_2=6.9,
\]
\[
\sum_{\text{ponderada}}=18+4(14.5)+2(6.9)=89.8,\qquad
c_s=\frac{h}{3}\sum_{\text{ponderada}}=\frac{89.8}{3}=29.93333333.
\]

Este desglose —tabla, \(n\), \(h\), sumatorias y resultado— se muestra en **Operaciones detalladas** y en el informe HTML exportado.

### D — Balance acumulado

La pendiente global \(k\) usa los puntos \((1,2.0),(3,5.3),(5,8.9),(7,12.2)\). También se calculan \(m\) con Simpson 1/3, \(n\) con bisección en \(1<n<3\), \(p\) con Runge-Kutta 4 desde \(y(0)=1\) y paso 0.1 hasta \(t=1.5\), y \(q\) con Newton hacia adelante en \(r=4.5\), pues la tabla de q tiene espaciamiento uniforme. Finalmente \(L(3)=k(3)^4+m(3)^3+n(3)^2+p(3)+q\) alimenta el parámetro \(\kappa\) de E.

### E — Dinámica del rendimiento

Se elige \(r_\sigma=1.2\), la menor raíz positiva de \(r^2-5.2r+4.8=0\), y se calcula \(\sigma=r_\sigma/10\). Se usa \(\kappa=L(3)\) del Documento D y \(\omega\) es el área de la tabla \(v(t)\) en \([0,3]\), calculada con Simpson 3/8. La ecuación se integra de 0 a 3 con Runge-Kutta 4 y paso 0.1; \(Y(0)=R(1.8)+x\).

**Supuesto declarado:** el expediente solo indica «acumulado 0–3» para \(\omega\) junto a la tabla \(v(t)\). Se calcula \(\omega=\int_0^3 v(t)\,dt\) con Simpson 3/8 sobre los cuatro puntos de la tabla, lo que da \(\omega=12\). Por tanto \(\kappa=L(3)\) y \(\omega\) son cantidades distintas y no se identifica \(\kappa/\omega\) con 1.

### F — Comportamiento térmico

Se calcula \(\chi\) por mínimos cuadrados, \(\psi=8\) como raíz mayor de \(r^2-11r+24=0\), y \(\xi=a+d+g\) del Documento A. Como el PDF no contiene una condición inicial legible para \(\theta\), se propone \(\theta(0)=25\,^\circ\mathrm C\), como temperatura ambiente o de reposo antes de iniciar la carga de trabajo. Se integra con Runge-Kutta 4 hasta \(t=3\) y paso 0.1; esta misma condición se usa en los análisis comparativos.

La ecuación original se conserva literalmente: \(\theta'=\chi\theta(\psi-\theta)-Y/\xi\). Se usa Runge-Kutta 4 y la trayectoria Y(t) del Documento E se interpola linealmente en cada etapa; \(\chi\), \(\psi\) y \(\xi\) son coeficientes constantes reconstruidos de sus datos de origen. La capacidad máxima del término logístico es \(\chi\psi^2/4\). Si \(Y/\xi\) supera ese máximo durante todo el intervalo, la derivada neta es negativa para cualquier \(\theta\); al caer \(\theta\) bajo cero, el término logístico también es negativo y la magnitud puede crecer hasta que RK4 deja de obtener valores finitos. El programa reporta esa explicación y el instante aproximado del fallo, sin inventar una temperatura final.

El signo negativo de \(Y/\xi\) hace que ese término reduzca la temperatura en la ecuación transcrita; no debe describirse como generación de calor que eleva la temperatura. La condición inicial de \(25\,^\circ\mathrm C\) es una propuesta justificada por el equipo, no un dato recuperado del PDF.

### G — Índice de estabilidad

Se resuelve \(a_sU^3-b_sU^2+c_sU-d_s=0\). Si hay varias raíces reales positivas, \(U^*\) es la menor. La integral de \(Y(t)\) en \([0,2]\) se calcula con Simpson 1/3 usando la trayectoria de E; \(\theta(1)\) se obtiene de la trayectoria de F. El índice solo se calcula si \(\theta(1)\) es finita. En esta instancia F diverge antes de \(t=1\), así que G y el índice no están disponibles; no se modifica F para forzar un resultado.

\[
S=\tfrac12\int_0^2Y(t)\,dt+\tfrac13\theta(1)+\ln(1+U^*).
\]

## Validación y sensibilidad

La pantalla **Análisis y gráficas** y el informe comparan Runge-Kutta 4 con \(h=0.1,0.05,0.025\), muestran las diferencias consecutivas cuando G está disponible, registran el estado/tiempo de fallo de F y repiten Newton–Raphson para C con tolerancias \(10^{-6},10^{-8},10^{-10}\). Las gráficas se construyen desde las trayectorias calculadas; F puede mostrarse parcialmente.

El menú **Operaciones detalladas** presenta sumatorias, iteraciones de Newton y bisección, matrices intermedias de Gauss–Jordan, el desglose de Simpson 1/3 para \(c_s\), resultados y supuestos. El botón **Exportar informe** guarda un HTML reproducible con tablas y gráficas SVG que incluye ese mismo desglose; se puede imprimir desde el navegador como PDF.

El HTML exportado incluye las matrículas de la instancia, resultados A–G, sensibilidad, validaciones numéricas, interpretación, supuestos y el contenido completo del archivo ejecutable. Imprimir ese HTML a PDF genera el informe técnico asociado a la instancia capturada.

## Limitaciones que deben reportarse

1. La condición inicial \(\theta(0)=25\,^\circ\mathrm C\) es propuesta por el equipo, no recuperada del archivo.
2. La selección de la menor raíz positiva para \(U^*\) es un criterio documentado para resolver la ambigüedad de selección.
3. La integración RK4 depende del paso \(h\); la comparación de pasos es un análisis empírico de sensibilidad, no una cota rigurosa del error global.
4. Conservar F sin normalización impide obtener G para esta instancia; cualquier resultado numérico de G requeriría justificar un cambio del modelo con el docente y el expediente.
