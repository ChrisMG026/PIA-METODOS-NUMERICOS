# Manual de usuario — PIA de Métodos Numéricos

## Requisitos

- Windows con Python 3.10 o posterior.
- Tkinter habilitado en la instalación de Python.
- No se instalan paquetes adicionales.

## Inicio

1. Abre una terminal en la carpeta del proyecto.
2. Ejecuta `python PIA.py`.
3. Captura cuatro matrículas distintas de exactamente siete dígitos y presiona **Calcular modelos**.
4. El programa las ordena automáticamente antes de asignar M1–M4. Conserva los ceros iniciales de la matrícula para mostrarla correctamente.

Si una matrícula no cumple el formato, aparece un mensaje para corregirla. Si algún cálculo numérico no puede completarse, la aplicación muestra el error en lugar de guardar un resultado incompleto.

## Navegación

- **Resultados** presenta una fila por cálculo y explica qué significa cada dato:
  - **M** es la matrícula escrita como un número.
  - **N** es ese número convertido a una escala pequeña para el modelo.
  - **P** es el promedio de los siete dígitos.
  - **S base** es la suma que da más peso a los dígitos ubicados al final.
- **VISTA GENERAL** contiene **Resultados base**, **Análisis y gráficas** y **Operaciones detalladas**; **DOCUMENTOS** contiene las vistas A–G.
- En **Operaciones detalladas** puedes seguir el proceso completo: iteraciones de Newton y bisección, matrices de Gauss-Jordan y el desglose de Simpson 1/3 para \(c_s\) (tabla, \(n\), \(h\), sumatorias y resultado). Ese mismo desglose aparece en el informe HTML.
- En cada documento A–G, la tabla muestra el resultado y una explicación breve. Las claves técnicas se conservan para relacionar las cifras con las ecuaciones del expediente.
- **Análisis y gráficas** compara la respuesta y la temperatura al usar distintos tamaños de paso de Runge-Kutta 4. Si el modelo térmico no llega al final, muestra hasta dónde pudo calcularse y deja el índice general como no disponible.
- Las tablas anchas se pueden recorrer con la barra de desplazamiento horizontal.
- **Exportar informe** permite guardar tablas, supuestos, análisis de sensibilidad, el desglose de Simpson 1/3 para \(c_s\) y gráficas en un archivo HTML.
- **Nueva captura** permite ingresar otra instancia del equipo.

El índice general **S** del Documento G es distinto de la suma ponderada **S base** de las matrículas. El índice general es un valor calculado, no una calificación: el expediente no define límites para decir que el sistema sea «estable» o «inestable», así que no debe interpretarse con esas etiquetas.

Para producir el entregable PDF, abre el HTML guardado en un navegador y selecciona **Imprimir → Guardar como PDF**. El HTML contiene los resultados de la instancia, las gráficas reproducibles, la validación y el código ejecutable completo.

## Supuestos que deben conservarse en el reporte

- \(\theta(0)=25\,^\circ\mathrm C\) es una condición inicial propuesta como temperatura ambiente/de reposo previa a la carga; declárala y justifícala en el reporte y úsala en todas las comparaciones.
- En la ecuación F transcrita, el término \(-Y/\xi\) reduce θ; no debe reportarse como generación de calor que eleva la temperatura.
- \(\kappa=L(3)\) se calcula con los coeficientes del Documento D; \(\omega\) se obtiene con Simpson 3/8 a partir de la tabla correspondiente.
- La ecuación F se conserva sin normalización. Para estas matrículas la trayectoria no alcanza un valor finito en todo el intervalo, por lo que G y el índice general aparecen como no disponibles. No se debe reemplazar ese estado por una cifra calculada con una ecuación modificada.
- Para el índice G se usa la menor raíz real positiva del modelo de saturación.

El informe técnico muestra el diagnóstico y el tiempo aproximado donde se detuvo F para que la limitación sea explícita.
