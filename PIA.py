import tkinter as tk
from tkinter import messagebox


def ln(x):
    """Calcula ln(x) mediante la serie de Mercator transformada."""
    # El logaritmo natural solo acepta argumentos positivos.
    if x <= 0:
        raise ValueError("El logaritmo natural solo está definido para x > 0.")
    # Esta sustitución permite usar la serie atanh(z) con mejor convergencia.
    z = (x - 1) / (x + 1)
    total = 0.0
    # El primer término de la serie es z/1.
    termino = z
    n = 1
    # Se acumulan términos hasta que su aporte sea menor que la tolerancia.
    while abs(termino) > 1e-12:
        total += termino / n
        # La serie usa potencias impares: z/1 + z^3/3 + z^5/5 + ...
        n += 2
        termino *= z * z
    # La relación ln(x) = 2 * (z + z^3/3 + z^5/5 + ...) completa el cálculo.
    return 2 * total


def exp(x):
    """Calcula exp(x) con la serie de Taylor de e^x."""
    total = 1.0
    termino = 1.0
    n = 1
    # Se suman los términos 1 + x + x^2/2! + x^3/3! + ...
    while abs(termino) > 1e-12:
        # Se obtiene el siguiente término usando el anterior, sin factoriales.
        termino *= x / n
        total += termino
        n += 1
    return total


def interpolacion_lagrange(x, nodos, valores):
    """Estima f(x) con el polinomio que pasa por los puntos dados."""
    resultado = 0.0
    for indice, nodo in enumerate(nodos):
        # Cada base de Lagrange comienza con el valor f(nodo).
        termino = valores[indice]
        for otro_indice, otro_nodo in enumerate(nodos):
            if indice != otro_indice:
                # Este cociente hace que la base valga 1 en su nodo y 0 en los demas.
                termino *= (x - otro_nodo) / (nodo - otro_nodo)
        # La interpolacion es la suma de todas las bases ponderadas.
        resultado += termino
    return resultado


def resolver_c_newton(P1, P2, S3, tolerancia=1e-10, max_iteraciones=100):
    """Resuelve la ecuacion no lineal de C con Newton-Raphson."""
    # La ecuacion de C contiene P1 + P2 varias veces, por eso se guarda una vez.
    suma_p = P1 + P2
    # Se elige una aproximacion inicial positiva y razonable para la iteracion.
    aproximacion = max(1.0, suma_p)
    for iteracion in range(1, max_iteraciones + 1):
        # f(u) = u^3 - (P1 + P2)u^2 + S3u - 5.
        funcion = (
            aproximacion ** 3
            - suma_p * aproximacion ** 2
            + S3 * aproximacion
            - 5
        )
        # f'(u) es la derivada necesaria para la formula de Newton-Raphson.
        derivada = 3 * aproximacion ** 2 - 2 * suma_p * aproximacion + S3
        # No se puede dividir entre una derivada nula o practicamente nula.
        if abs(derivada) < 1e-14:
            raise ValueError("La derivada se hizo cero durante Newton-Raphson.")
        # Formula: u_nuevo = u_actual - f(u_actual) / f'(u_actual).
        nueva_aproximacion = aproximacion - funcion / derivada
        # Si el cambio ya es pequeño, se considera que se encontro la raiz.
        if abs(nueva_aproximacion - aproximacion) < tolerancia:
            return nueva_aproximacion, iteracion
        aproximacion = nueva_aproximacion
    # Evita un ciclo indefinido si el metodo no converge.
    raise ValueError("Newton-Raphson no convergio en el limite de iteraciones.")


def calcular_coeficientes_documento_b():
    """Calcula Ar, Br, Cr y Dr para el modelo de respuesta R(n)."""
    def funcion(w):
        return w ** 3 - 2.1 * w ** 2 + 0.9 * w - 0.08

    def derivada(w):
        return 3 * w ** 2 - 4.2 * w + 0.9

    raiz = 0.1
    for _ in range(100):
        valor_derivada = derivada(raiz)
        if abs(valor_derivada) < 1e-14:
            raise ValueError("No se pudo calcular Ar: derivada nula.")
        siguiente = raiz - funcion(raiz) / valor_derivada
        if abs(siguiente - raiz) < 1e-12:
            raiz = siguiente
            break
        raiz = siguiente
    else:
        raise ValueError("No hubo convergencia al calcular Ar.")

    Ar = raiz
    Br = 2 * ln(3) + 13 / 3
    nodos = (1, 2, 4, 6)
    valores = (4.2, 6.1, 9.8, 13.9)
    media_n = sum(nodos) / len(nodos)
    media_h = sum(valores) / len(valores)
    Cr = sum(
        (nodo - media_n) * (valor - media_h)
        for nodo, valor in zip(nodos, valores)
    ) / sum((nodo - media_n) ** 2 for nodo in nodos)
    Dr = media_h - Cr * media_n
    return {"Ar": Ar, "Br": Br, "Cr": Cr, "Dr": Dr}


def resolver_gauss_jordan(matriz, tolerancia=1e-12, pasos=None):
    """Resuelve un sistema cuadrado mediante Gauss-Jordan con pivoteo parcial y opcionalmente guarda pasos."""
    sistema = [fila[:] for fila in matriz]
    cantidad = len(sistema)
    for columna in range(cantidad):
        fila_pivote = max(
            range(columna, cantidad),
            key=lambda fila: abs(sistema[fila][columna]),
        )
        sistema[columna], sistema[fila_pivote] = sistema[fila_pivote], sistema[columna]
        pivote = sistema[columna][columna]
        if abs(pivote) < tolerancia:
            raise ValueError("El sistema de Documento A no tiene solución única.")
        sistema[columna] = [valor / pivote for valor in sistema[columna]]
        for fila in range(cantidad):
            if fila != columna:
                factor = sistema[fila][columna]
                sistema[fila] = [
                    actual - factor * pivote_actual
                    for actual, pivote_actual in zip(sistema[fila], sistema[columna])
                ]
        if pasos is not None:
            pasos.append([fila[:] for fila in sistema])
    return [sistema[fila][-1] for fila in range(cantidad)]


def calcular_documento_c(variables_a):
    """Calcula a_s y b_s con los coeficientes de Documento B y A."""
    coeficientes = calcular_coeficientes_documento_b()
    Ar, Br, Cr, Dr = (coeficientes[nombre] for nombre in ("Ar", "Br", "Cr", "Dr"))
    a_s = Ar * 2 ** 3 + Br * 2 ** 2 + Cr * 2 + Dr

    matriz = [
        [variables_a[fila]["valor"] for fila in "abc"]
        + [variables_a["A"]["valor"]],
        [variables_a[fila]["valor"] for fila in "def"]
        + [variables_a["B"]["valor"]],
        [variables_a[fila]["valor"] for fila in "ghi"]
        + [variables_a["C"]["valor"]],
    ]
    x, y, z = resolver_gauss_jordan(matriz)
    b_s = x + y + z
    return {
        "Ar": Ar,
        "Br": Br,
        "Cr": Cr,
        "Dr": Dr,
        "a_s": a_s,
        "x": x,
        "y": y,
        "z": z,
        "b_s": b_s,
    }


def calcular_resultados(matriculas):
    """Calcula M, N, P y S a partir de las matriculas capturadas."""
    M = {
        clave: int("".join(str(digito) for digito in digitos))
        for clave, digitos in matriculas.items()
    }
    N = {clave: valor / (10 ** 7) for clave, valor in M.items()}
    P = {
        clave: (1 / 7) * sum(digitos)
        for clave, digitos in matriculas.items()
    }
    S = {
        clave: sum(pos * digito for pos, digito in enumerate(digitos, start=1))
        for clave, digitos in matriculas.items()
    }
    return M, N, P, S


class Aplicacion(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Métodos numéricos - Matrículas")
        self.geometry("900x700")
        self.minsize(720, 600)
        self.configure(bg="#eef2f7")
        self.matriculas = {}
        self.resultados = None
        self.entradas = []
        self.pantalla_actual = "captura"
        self.mostrar_captura()

    def limpiar_pantalla(self):
        for widget in self.winfo_children():
            widget.destroy()

    def crear_canvas(self, titulo):
        self.limpiar_pantalla()
        canvas = tk.Canvas(self, bg="#eef2f7", highlightthickness=0)
        canvas.pack(fill="both", expand=True)
        canvas.create_text(
            45, 40, anchor="w", text=titulo,
            font=("Segoe UI", 24, "bold"), fill="#172033"
        )
        return canvas

    def mostrar_captura(self):
        self.pantalla_actual = "captura"
        canvas = self.crear_canvas("Captura de matrículas")
        canvas.create_text(
            45, 82, anchor="w",
            text="Escribe cuatro matrículas distintas de 7 dígitos.",
            font=("Segoe UI", 11), fill="#526074"
        )
        panel = tk.Frame(canvas, bg="white", padx=28, pady=20)
        canvas.create_window(45, 115, anchor="nw", window=panel, width=500)
        self.entradas = []
        for indice in range(1, 5):
            tk.Label(
                panel, text=f"M{indice}", bg="white", fg="#172033",
                font=("Segoe UI", 11, "bold")
            ).grid(row=indice - 1, column=0, padx=(0, 18), pady=8)
            entrada = tk.Entry(
                panel, width=18, font=("Consolas", 13),
                justify="center", relief="solid", bd=1
            )
            entrada.grid(row=indice - 1, column=1, pady=8)
            self.entradas.append(entrada)
        boton = tk.Button(
            panel, text="Calcular resultados", command=self.procesar,
            bg="#216e9b", fg="white", activebackground="#185574",
            activeforeground="white", relief="flat", padx=16, pady=9,
            font=("Segoe UI", 10, "bold"), cursor="hand2"
        )
        boton.grid(row=4, column=0, columnspan=2, pady=(18, 0))

    def procesar(self):
        valores = [entrada.get().strip() for entrada in self.entradas]
        if any(not valor.isdigit() or len(valor) != 7 for valor in valores):
            messagebox.showerror(
                "Dato inválido",
                "Cada matrícula debe contener exactamente 7 dígitos numéricos."
            )
            return
        if len(set(valores)) != 4:
            messagebox.showerror("Dato inválido", "Las matrículas no pueden repetirse.")
            return
        self.matriculas = {
            f"M{indice}": [int(digito) for digito in valor]
            for indice, valor in enumerate(valores, start=1)
        }
        self.resultados = calcular_resultados(self.matriculas)
        self.mostrar_resultados()

    def calcular_documento_a(self):
        if not self.resultados:
            return {}
        _, N, P, S = self.resultados
        a = P["M1"] + N["M2"]
        b = S["M2"] ** 0.5
        c = (P["M3"] ** 2) / P["M4"]
        d = ln(1 + S["M1"])
        e = (P["M2"] ** 2) + N["M3"]
        f = S["M4"] / 10
        g = exp(N["M1"])
        h = S["M3"] / (P["M1"] + 1)
        i = (P["M4"] + S["M2"]) ** 0.5
        A = (2 / 3) * (P["M1"] ** 3) + (S["M2"] / 2) * (P["M1"] ** 2) + P["M1"]
        nodos_b = (P["M2"], P["M3"], P["M4"])
        valores_b = (4, 9, 16)
        B = interpolacion_lagrange(N["M2"], nodos_b, valores_b)
        C, iteraciones_c = resolver_c_newton(P["M1"], P["M2"], S["M3"])
        D = interpolacion_lagrange(2.5, (1, 2, 3, 4), (2.8, 4.6, 7.9, 12.5))
        variables = {
            "a": {"formula": "P1 + N2", "valor": a},
            "b": {"formula": "sqrt(S2)", "valor": b},
            "c": {"formula": "P3^2 / P4", "valor": c},
            "d": {"formula": "ln(1 + S1)", "valor": d},
            "e": {"formula": "P2^2 + N3", "valor": e},
            "f": {"formula": "S4 / 10", "valor": f},
            "g": {"formula": "exp(N1)", "valor": g},
            "h": {"formula": "S3 / (P1 + 1)", "valor": h},
            "i": {"formula": "sqrt(P4 + S2)", "valor": i},
            "A": {"formula": "(2/3)P1^3 + (S2/2)P1^2 + P1", "valor": A},
            "B": {"formula": "x=N2; (P2,4), (P3,9), (P4,16)", "valor": B},
            "C": {"formula": "u^3 - (P1+P2)u^2 + S3u - 5 = 0", "valor": C},
            "D": {"formula": "x=2.5; tabla de interpolacion", "valor": D},
        }
        return variables

    def agregar_navegacion(self, canvas):
        barra = tk.Frame(self, bg="#172033", height=52)
        barra.place(relx=0, rely=1, relwidth=1, anchor="sw")
        tk.Button(
            barra, text="Resultados", command=self.mostrar_resultados,
            bg="#172033", fg="white", relief="flat", padx=18,
            font=("Segoe UI", 10, "bold")
        ).pack(side="left", padx=18, pady=10)
        tk.Button(
            barra, text="Documento A", command=self.mostrar_documento_a,
            bg="#172033", fg="white", relief="flat", padx=18,
            font=("Segoe UI", 10, "bold")
        ).pack(side="left", pady=10)
        tk.Button(
            barra, text="Documento B", command=self.mostrar_documento_b,
            bg="#172033", fg="white", relief="flat", padx=18,
            font=("Segoe UI", 10, "bold")
        ).pack(side="left", pady=10)
        tk.Button(
            barra, text="Documento C", command=self.mostrar_documento_c,
            bg="#172033", fg="white", relief="flat", padx=18,
            font=("Segoe UI", 10, "bold")
        ).pack(side="left", pady=10)
        tk.Button(
            barra, text="Operaciones e iteraciones", command=self.mostrar_operaciones,
            bg="#172033", fg="white", relief="flat", padx=18,
            font=("Segoe UI", 10, "bold")
        ).pack(side="left", pady=10)
        tk.Button(
            barra, text="Nueva captura", command=self.mostrar_captura,
            bg="#172033", fg="#9bd4f3", relief="flat", padx=18,
            font=("Segoe UI", 10, "bold")
        ).pack(side="right", padx=18, pady=10)

    def mostrar_resultados(self):
        if not self.resultados:
            self.mostrar_captura()
            return
        self.pantalla_actual = "resultados"
        canvas = self.crear_canvas("Resultados")
        M, N, P, S = self.resultados
        encabezados = ("Matrícula", "M", "N", "P", "S")
        columnas = (120, 260, 400, 540, 680)
        for x, encabezado in zip(columnas, encabezados):
            canvas.create_text(x, 105, text=encabezado, font=("Segoe UI", 11, "bold"), fill="#216e9b")
        canvas.create_line(55, 122, 830, 122, fill="#b8c4d3")
        for fila, clave in enumerate(("M1", "M2", "M3", "M4")):
            y = 160 + fila * 62
            matricula = "".join(str(digito) for digito in self.matriculas[clave])
            valores = (matricula, str(M[clave]), f"{N[clave]:.7f}", f"{P[clave]:.7f}", str(S[clave]))
            for x, valor in zip(columnas, valores):
                canvas.create_text(x, y, text=valor, font=("Consolas", 12), fill="#172033")
            canvas.create_line(55, y + 25, 830, y + 25, fill="#d7dee8")
        self.agregar_navegacion(canvas)

    def mostrar_documento_a(self):
        if not self.resultados:
            self.mostrar_captura()
            return
        self.pantalla_actual = "documento_a"
        canvas = self.crear_canvas("Documento A: Estado Operativo")
        variables = self.calcular_documento_a()

        canvas.create_text(100, 90, anchor="w", text="Variable", font=("Segoe UI", 11, "bold"), fill="#216e9b")
        canvas.create_text(420, 90, anchor="w", text="Valor", font=("Segoe UI", 11, "bold"), fill="#216e9b")
        canvas.create_line(45, 105, 840, 105, fill="#b8c4d3")

        letras = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "A", "B", "C", "D"]
        for indice, letra in enumerate(letras):
            y = 130 + indice * 38
            canvas.create_text(100, y, anchor="w", text=letra, font=("Consolas", 12), fill="#172033")
            canvas.create_text(420, y, anchor="w", text=f"{variables[letra]['valor']:.7f}", font=("Consolas", 12), fill="#172033")
            canvas.create_line(45, y + 18, 840, y + 18, fill="#d7dee8")
        self.agregar_navegacion(canvas)

    def mostrar_documento_b(self):
        self.pantalla_actual = "documento_b"
        canvas = self.crear_canvas("Documento B: Respuesta bajo carga")
        coeficientes = calcular_coeficientes_documento_b()

        canvas.create_text(150, 90, anchor="w", text="Variable", font=("Segoe UI", 11, "bold"), fill="#216e9b")
        canvas.create_text(450, 90, anchor="w", text="Valor", font=("Segoe UI", 11, "bold"), fill="#216e9b")
        canvas.create_line(45, 105, 840, 105, fill="#b8c4d3")
        for indice, nombre in enumerate(("Ar", "Br", "Cr", "Dr")):
            y = 145 + indice * 52
            canvas.create_text(150, y, anchor="w", text=nombre, font=("Consolas", 12), fill="#172033")
            canvas.create_text(450, y, anchor="w", text=f"{coeficientes[nombre]:.10f}", font=("Consolas", 12), fill="#172033")
            canvas.create_line(45, y + 20, 840, y + 20, fill="#d7dee8")
        self.agregar_navegacion(canvas)

    def mostrar_documento_c(self):
        if not self.resultados:
            self.mostrar_captura()
            return
        self.pantalla_actual = "documento_c"
        canvas = self.crear_canvas("Documento C: Saturación")
        variables_a = self.calcular_documento_a()
        canvas.create_text(
            45, 78, anchor="w",
            text="Coeficientes disponibles del modelo: a_s y b_s",
            font=("Segoe UI", 13, "bold"), fill="#216e9b"
        )

        def calcular():
            try:
                valores = calcular_documento_c(variables_a)
            except ValueError as error:
                messagebox.showerror("No se pudo calcular Documento C", str(error))
                return
            return valores

        try:
            valores = calcular()
            filas = (("a_s", valores["a_s"]), ("b_s", valores["b_s"]))
        except (TypeError, KeyError):
            filas = (("a_s", "No disponible"), ("b_s", "No disponible"))

        canvas.create_text(150, 100, anchor="w", text="Variable", font=("Segoe UI", 11, "bold"), fill="#216e9b")
        canvas.create_text(450, 100, anchor="w", text="Valor", font=("Segoe UI", 11, "bold"), fill="#216e9b")
        canvas.create_line(45, 115, 840, 115, fill="#b8c4d3")
        for indice, (nombre, valor) in enumerate(filas):
            y = 155 + indice * 52
            texto_valor = f"{valor:.7f}" if isinstance(valor, (int, float)) else valor
            canvas.create_text(150, y, anchor="w", text=nombre, font=("Consolas", 12), fill="#172033")
            canvas.create_text(450, y, anchor="w", text=texto_valor, font=("Consolas", 12), fill="#172033")
            canvas.create_line(45, y + 20, 840, y + 20, fill="#d7dee8")
        self.agregar_navegacion(canvas)

    def mostrar_operaciones(self):
        if not self.resultados:
            self.mostrar_captura()
            return
        self.pantalla_actual = "operaciones"
        canvas = self.crear_canvas("Operaciones e iteraciones")
        M, N, P, S = self.resultados
        variables_a = self.calcular_documento_a()
        variables_c = calcular_documento_c(variables_a)
        _, iteraciones_c = resolver_c_newton(P["M1"], P["M2"], S["M3"])
        residuo_c = (
            variables_a["C"]["valor"] ** 3
            - (P["M1"] + P["M2"]) * variables_a["C"]["valor"] ** 2
            + S["M3"] * variables_a["C"]["valor"]
            - 5
        )
        matriz_c = [
            [variables_a[fila]["valor"] for fila in "abc"] + [variables_a["A"]["valor"]],
            [variables_a[fila]["valor"] for fila in "def"] + [variables_a["B"]["valor"]],
            [variables_a[fila]["valor"] for fila in "ghi"] + [variables_a["C"]["valor"]],
        ]
        pasos_gauss = []
        resolver_gauss_jordan(matriz_c, pasos=pasos_gauss)
        lineas = [
            "DATOS DE MATRÍCULAS",
            "Validación: cuatro matrículas distintas de siete dígitos.",
            "",
        ]
        for clave in ("M1", "M2", "M3", "M4"):
            digitos = self.matriculas[clave]
            suma_digitos = " + ".join(str(digito) for digito in digitos)
            terminos_s = " + ".join(f"{pos}({digito})" for pos, digito in enumerate(digitos, start=1))
            lineas.extend([
                f"{clave}: M = int({''.join(str(digito) for digito in digitos)}) = {M[clave]}",
                f"  N = M / 10^7 = {M[clave]} / 10000000 = {N[clave]:.7f}",
                f"  P = (1/7)({suma_digitos}) = {P[clave]:.7f}",
                f"  S = {terminos_s} = {S[clave]}",
                "",
            ])

        lineas.extend(["DOCUMENTO A: ESTADO OPERATIVO", ""])
        lineas.extend([
            f"a = P1 + N2 = {P['M1']:.7f} + {N['M2']:.7f} = {variables_a['a']['valor']:.7f}",
            f"b = sqrt(S2) = sqrt({S['M2']}) = {variables_a['b']['valor']:.7f}",
            f"c = P3^2 / P4 = {P['M3']:.7f}^2 / {P['M4']:.7f} = {variables_a['c']['valor']:.7f}",
            f"d = ln(1 + S1) = ln(1 + {S['M1']}) = {variables_a['d']['valor']:.7f}",
            f"e = P2^2 + N3 = {P['M2']:.7f}^2 + {N['M3']:.7f} = {variables_a['e']['valor']:.7f}",
            f"f = S4 / 10 = {S['M4']} / 10 = {variables_a['f']['valor']:.7f}",
            f"g = exp(N1) = exp({N['M1']:.7f}) = {variables_a['g']['valor']:.7f}",
            f"h = S3 / (P1 + 1) = {S['M3']} / ({P['M1']:.7f} + 1) = {variables_a['h']['valor']:.7f}",
            f"i = sqrt(P4 + S2) = sqrt({P['M4']:.7f} + {S['M2']}) = {variables_a['i']['valor']:.7f}",
            f"A = (2/3)P1^3 + (S2/2)P1^2 + P1 = {variables_a['A']['valor']:.7f}",
            "B = interpolación de Lagrange con x=N2 y los nodos P2, P3, P4:",
            f"  valores tabulados (4, 9, 16); B = {variables_a['B']['valor']:.7f}",
            "C = raíz por Newton-Raphson de F(u)=u^3-(P1+P2)u^2+S3u-5=0.",
            f"  F'(u)=3u^2-2(P1+P2)u+S3; C={variables_a['C']['valor']:.7f}",
            f"  Iteraciones={iteraciones_c}; residuo F(C)={residuo_c:.3e}",
            "D = interpolación de Lagrange evaluada en x=2.5 con (1,2.8), (2,4.6), (3,7.9), (4,12.5):",
            f"  D = {variables_a['D']['valor']:.7f}",
            "",
        ])

        coeficientes_b = calcular_coeficientes_documento_b()
        lineas.extend([
            "DOCUMENTO B: RESPUESTA BAJO CARGA",
            "",
            "Ar: Newton-Raphson en w^3 - 2.1w^2 + 0.9w - 0.08 = 0.",
            f"  Ar = {coeficientes_b['Ar']:.10f}",
            "Br: integral de 1 a 3 de (2/x + 0.5x^2) dx = 2ln(3) + 13/3.",
            f"  Br = {coeficientes_b['Br']:.10f}",
            "Cr y Dr: tendencia lineal global con (1,4.2), (2,6.1), (4,9.8), (6,13.9).",
            f"  Cr = {coeficientes_b['Cr']:.10f}; Dr = {coeficientes_b['Dr']:.10f}",
            "",
            "DOCUMENTO C: SATURACIÓN (COEFICIENTES DISPONIBLES)",
            "",
            f"a_s = R(2) = Ar(2^3) + Br(2^2) + Cr(2) + Dr = {variables_c['a_s']:.7f}",
            "Sistema para Gauss-Jordan: ax+by+cz=A; dx+ey+fz=B; gx+hy+iz=C.",
            "Matriz aumentada:",
            f"  [{variables_a['a']['valor']:.5f}, {variables_a['b']['valor']:.5f}, {variables_a['c']['valor']:.5f} | {variables_a['A']['valor']:.5f}]",
            f"  [{variables_a['d']['valor']:.5f}, {variables_a['e']['valor']:.5f}, {variables_a['f']['valor']:.5f} | {variables_a['B']['valor']:.5f}]",
            f"  [{variables_a['g']['valor']:.5f}, {variables_a['h']['valor']:.5f}, {variables_a['i']['valor']:.5f} | {variables_a['C']['valor']:.5f}]",
            "Matrices tras cada pivote de Gauss-Jordan:",
            *[
                f"  Paso {indice}: " + " ; ".join(
                    "[" + ", ".join(f"{valor:.5f}" for valor in fila) + "]"
                    for fila in paso
                )
                for indice, paso in enumerate(pasos_gauss, start=1)
            ],
            f"Solución: x={variables_c['x']:.7f}, y={variables_c['y']:.7f}, z={variables_c['z']:.7f}.",
            f"b_s = x + y + z = {variables_c['b_s']:.7f}",
            "",
        ])

        contenedor = tk.Frame(self, bg="white")
        contenedor.place(relx=0.05, rely=0.13, relwidth=0.9, relheight=0.78)
        texto = tk.Text(
            contenedor, bg="white", fg="#172033", font=("Consolas", 10),
            relief="flat", padx=18, pady=14, wrap="word"
        )
        texto.insert("1.0", "\n".join(lineas))
        texto.configure(state="disabled")
        scrollbar = tk.Scrollbar(contenedor, command=texto.yview)
        texto.configure(yscrollcommand=scrollbar.set)
        texto.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        self.agregar_navegacion(canvas)


if __name__ == "__main__":
    Aplicacion().mainloop()