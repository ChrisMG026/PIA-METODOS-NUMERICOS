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
            barra, text="Documento A: Estado Operativo", command=self.mostrar_documento_a,
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

        canvas.create_text(50, 90, anchor="w", text="Variable", font=("Segoe UI", 11, "bold"), fill="#216e9b")
        canvas.create_text(230, 90, anchor="w", text="Fórmula", font=("Segoe UI", 11, "bold"), fill="#216e9b")
        canvas.create_text(500, 90, anchor="w", text="Valor", font=("Segoe UI", 11, "bold"), fill="#216e9b")
        canvas.create_line(45, 105, 840, 105, fill="#b8c4d3")

        letras = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "A", "B", "C"]
        for indice, letra in enumerate(letras):
            y = 140 + indice * 42
            detalle = variables.get(letra, {"formula": "Pendiente", "valor": "—"})
            canvas.create_text(50, y, anchor="w", text=letra, font=("Consolas", 12), fill="#172033")
            canvas.create_text(230, y, anchor="w", text=detalle["formula"], font=("Consolas", 12), fill="#172033")
            canvas.create_text(500, y, anchor="w", text=str(detalle["valor"]), font=("Consolas", 12), fill="#172033")
            canvas.create_line(45, y + 18, 840, y + 18, fill="#d7dee8")

        canvas.create_text(
            50, 590,
            anchor="w",
            text="",
            font=("Segoe UI", 10),
            fill="#526074"
        )
        self.agregar_navegacion(canvas)

    def mostrar_operaciones(self):
        if not self.resultados:
            self.mostrar_captura()
            return
        self.pantalla_actual = "operaciones"
        canvas = self.crear_canvas("Operaciones e iteraciones")
        M, N, P, S = self.resultados
        lineas = [
            "Iteración de captura: se validaron 4 matrículas, todas distintas y de 7 dígitos.",
            "",
        ]
        for clave in ("M1", "M2", "M3", "M4"):
            digitos = self.matriculas[clave]
            suma_digitos = " + ".join(str(digito) for digito in digitos)
            terminos_s = " + ".join(f"{pos}({digito})" for pos, digito in enumerate(digitos, start=1))
            lineas.extend([
                f"{clave}: M = int({''.join(str(digito) for digito in digitos)}) = {M[clave]}",
                f"     N = M / 10^7 = {M[clave]} / 10000000 = {N[clave]:.7f}",
                f"     P = (1/7) * ({suma_digitos}) = {P[clave]:.7f}",
                f"     S = {terminos_s} = {S[clave]}",
                "",
            ])

        a = P["M1"] + N["M2"]
        b = S["M2"] ** 0.5
        c = (P["M3"] ** 2) / P["M4"]
        d = ln(1 + S["M1"])
        e = (P["M2"] ** 2) + N["M3"]
        f = S["M4"] / 10
        g = exp(N["M1"])
        h = S["M3"] / (P["M1"] + 1)
        i = (P["M4"] + S["M2"]) ** 0.5
        B = interpolacion_lagrange(
            N["M2"],
            (P["M2"], P["M3"], P["M4"]),
            (4, 9, 16),
        )
        C, iteraciones_c = resolver_c_newton(P["M1"], P["M2"], S["M3"])
        residuo_c = (
            C ** 3
            - (P["M1"] + P["M2"]) * C ** 2
            + S["M3"] * C
            - 5
        )
        lineas.extend([
            "Documento A: Estado Operativo",
            "",
            f"a = P1 + N2 = {P['M1']:.7f} + {N['M2']:.7f} = {a:.7f}",
            f"b = sqrt(S2) = sqrt({S['M2']}) = {b:.7f}",
            f"c = P3^2 / P4 = ({P['M3']:.7f})^2 / {P['M4']:.7f} = {c:.7f}",
            f"d = ln(1 + S1) = ln(1 + {S['M1']}) = {d:.7f}",
            f"e = P2^2 + N3 = ({P['M2']:.7f})^2 + {N['M3']:.7f} = {e:.7f}",
            f"f = S4 / 10 = {S['M4']} / 10 = {f:.7f}",
            f"g = exp(N1) = exp({N['M1']:.7f}) = {g:.7f}",
            f"h = S3 / (P1 + 1) = {S['M3']} / ({P['M1']:.7f} + 1) = {h:.7f}",
            f"i = sqrt(P4 + S2) = sqrt({P['M4']:.7f} + {S['M2']}) = {i:.7f}",
            "B: interpolacion de Lagrange con x = N2 y los puntos (P2,4), (P3,9), (P4,16)",
            f"B = {B:.7f}",
            "C: Newton-Raphson para u^3 - (P1 + P2)u^2 + S3u - 5 = 0",
            f"C = {C:.7f} despues de {iteraciones_c} iteraciones",
            f"Residuo de C = {residuo_c:.3e}",
            "",
        ])

        texto = tk.Text(
            self, bg="white", fg="#172033", font=("Consolas", 11),
            relief="flat", padx=22, pady=18, wrap="none"
        )
        texto.insert("1.0", "\n".join(lineas))
        texto.configure(state="disabled")
        texto.place(relx=0.05, rely=0.17, relwidth=0.9, relheight=0.68)
        self.agregar_navegacion(canvas)


if __name__ == "__main__":
    Aplicacion().mainloop()