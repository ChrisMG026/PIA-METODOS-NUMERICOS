# Manual de usuario — PIA de Métodos Numéricos

## Requisitos

- Windows con Python 3.10 o posterior.
- Tkinter habilitado en la instalación de Python.
- No se instalan paquetes adicionales.

## Inicio

1. Abre una terminal en la carpeta del proyecto.
2. Ejecuta `python PIA.py`.
3. Captura cuatro matrículas distintas de exactamente siete dígitos y presiona **Calcular Modelos Numéricos**.
4. El programa las ordena automáticamente antes de asignar M1–M4. Conserva los ceros iniciales de la matrícula para mostrarla correctamente.

Si una matrícula no cumple el formato, aparece un mensaje para corregirla. Si algún cálculo numérico no puede completarse, la aplicación muestra el error en lugar de guardar un resultado incompleto.

## Navegación

- **Resultados** presenta una fila por cálculo y explica qué significa cada dato:
  - **M** es la matrícula escrita como un número.
  - **N** es ese número convertido a una escala pequeña para el modelo.
  - **P** es el promedio de los siete dígitos.
  - **S base** es la suma que da más peso a los dígitos ubicados al final.
- **Documentos** contiene las vistas A–G y **Operaciones detalladas**.
- En cada documento A–G, la tabla muestra el resultado y una explicación breve. Las claves técnicas se conservan para relacionar las cifras con las ecuaciones del expediente.
- **Análisis y gráficas** compara la respuesta, la temperatura y el índice general al usar distintos tamaños de paso. En las gráficas, el eje horizontal representa el tiempo.
- Las tablas anchas se pueden recorrer con la barra de desplazamiento horizontal.
- **Exportar informe** permite guardar tablas, supuestos, análisis de sensibilidad y gráficas en un archivo HTML.
- **Nueva captura** permite ingresar otra instancia del equipo.

El índice general **S** del Documento G es distinto de la suma ponderada **S base** de las matrículas. El índice general es un valor calculado, no una calificación: el expediente no define límites para decir que el sistema sea «estable» o «inestable», así que no debe interpretarse con esas etiquetas.

Para producir el entregable PDF, abre el HTML guardado en un navegador y selecciona **Imprimir → Guardar como PDF**. El HTML contiene los resultados de la instancia, las gráficas reproducibles, la validación y el código ejecutable completo.

## Supuestos que deben conservarse en el reporte

- \(\theta(0)=25\,^\circ\mathrm C\) se propone como temperatura ambiente/de reposo previa al inicio de la carga.
- En la ecuación F transcrita, el término \(-Y/\xi\) reduce θ; no debe reportarse como generación de calor que eleva la temperatura.
- \(\omega=\kappa\) se usa porque el PDF no proporciona un valor separado para \(\omega\).
- La ecuación térmica del expediente diverge con datos típicos. Para obtener una trayectoria finita, el programa aplica el factor \(\lambda\) descrito en el manual técnico, limitando la fuente a 80 % de la capacidad logística. Esta corrección no aparece en el PDF y debe mencionarse como una decisión del equipo.
- Para el índice G se usa la menor raíz real positiva del modelo de saturación.

El usuario debe incluir estas decisiones y sus razones en el manual técnico y no presentar la normalización térmica como parte original del expediente.
