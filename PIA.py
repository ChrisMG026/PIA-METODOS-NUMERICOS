import tkinter as tk
from tkinter import messagebox


def ln(x):
    """Calcula ln(x) mediante la serie de Mercator transformada."""
    if x <= 0:
        raise ValueError("El logaritmo natural solo está definido para x > 0.")
    z = (x - 1) / (x + 1)
    total = 0.0
    termino = z
    n = 1
    while abs(termino) > 1e-12:
        total += termino / n
        n += 2
        termino *= z * z
    return 2 * total


def exp(x):
    """Calcula exp(x) con la serie de Taylor de e^x."""
    total = 1.0
    termino = 1.0
    n = 1
    while abs(termino) > 1e-12:
        termino *= x / n
        total += termino
        n += 1
    return total


def interpolacion_lagrange(x, nodos, valores):
    """Estima f(x) con el polinomio que pasa por los puntos dados."""
    resultado = 0.0
    for indice, nodo in enumerate(nodos):
        termino = valores[indice]
        for otro_indice, otro_nodo in enumerate(nodos):
            if indice != otro_indice:
                termino *= (x - otro_nodo) / (nodo - otro_nodo)
        resultado += termino
    return resultado


def resolver_c_newton(P1, P2, S3, tolerancia=1e-10, max_iteraciones=100):
    """Resuelve la ecuacion no lineal de C con Newton-Raphson registrando iteraciones."""
    suma_p = P1 + P2
    aproximacion = max(1.0, suma_p)
    historial = []
    for iteracion in range(1, max_iteraciones + 1):
        funcion = (
            aproximacion ** 3
            - suma_p * aproximacion ** 2
            + S3 * aproximacion
            - 5
        )
        derivada = 3 * aproximacion ** 2 - 2 * suma_p * aproximacion + S3
        if abs(derivada) < 1e-14:
            raise ValueError("La derivada se hizo cero durante Newton-Raphson.")
        nueva_aproximacion = aproximacion - funcion / derivada
        historial.append({
            "iter": iteracion,
            "x": aproximacion,
            "fx": funcion,
            "dfx": derivada,
            "x_nueva": nueva_aproximacion
        })
        if abs(nueva_aproximacion - aproximacion) < tolerancia:
            return nueva_aproximacion, iteracion, historial
        aproximacion = nueva_aproximacion
    raise ValueError("Newton-Raphson no convergio en el limite de iteraciones.")


def calcular_coeficientes_documento_b():
    """Calcula Ar, Br, Cr y Dr para el modelo de respuesta R(n)."""
    def funcion(w):
        return w ** 3 - 2.1 * w ** 2 + 0.9 * w - 0.08

    def derivada(w):
        return 3 * w ** 2 - 4.2 * w + 0.9

    raiz = 0.1
    historial_ar = []
    for i in range(1, 100):
        val_f = funcion(raiz)
        val_df = derivada(raiz)
        if abs(val_df) < 1e-14:
            raise ValueError("No se pudo calcular Ar: derivada nula.")
        siguiente = raiz - val_f / val_df
        historial_ar.append({"iter": i, "x": raiz, "fx": val_f, "dfx": val_df, "x_nueva": siguiente})
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
    return {"Ar": Ar, "Br": Br, "Cr": Cr, "Dr": Dr, "historial_ar": historial_ar}


def resolver_gauss_jordan(matriz, tolerancia=1e-12):
    """Resuelve un sistema cuadrado mediante Gauss-Jordan con registro detallado de pasos."""
    sistema = [fila[:] for fila in matriz]
    cantidad = len(sistema)
    pasos = []
    pasos.append({"etapa": "Matriz Aumentada Inicial", "matriz": [fila[:] for fila in sistema]})

    for columna in range(cantidad):
        fila_pivote = max(
            range(columna, cantidad),
            key=lambda fila: abs(sistema[fila][columna]),
        )
        if fila_pivote != columna:
            sistema[columna], sistema[fila_pivote] = sistema[fila_pivote], sistema[columna]
            pasos.append({
                "etapa": f"Intercambio de fila {columna+1} con fila {fila_pivote+1}",
                "matriz": [fila[:] for fila in sistema]
            })

        pivote = sistema[columna][columna]
        if abs(pivote) < tolerancia:
            raise ValueError("El sistema no tiene solución única.")

        sistema[columna] = [valor / pivote for valor in sistema[columna]]
        pasos.append({
            "etapa": f"Normalizar fila {columna+1} dividiendo entre el pivote {pivote:.5f}",
            "matriz": [fila[:] for fila in sistema]
        })

        for fila in range(cantidad):
            if fila != columna:
                factor = sistema[fila][columna]
                sistema[fila] = [
                    actual - factor * pivote_actual
                    for actual, pivote_actual in zip(sistema[fila], sistema[columna])
                ]
                pasos.append({
                    "etapa": f"Eliminar elemento en fila {fila+1}, columna {columna+1} (factor {factor:.5f})",
                    "matriz": [fila[:] for fila in sistema]
                })

    return [sistema[fila][-1] for fila in range(cantidad)], pasos


def calcular_documento_c(variables_a):
    """Calcula a_s, b_s, c_s y d_s para el Documento C."""
    coeficientes = calcular_coeficientes_documento_b()
    Ar, Br, Cr, Dr = (coeficientes[nombre] for nombre in ("Ar", "Br", "Cr", "Dr"))
    a_s = Ar * (2 ** 3) + Br * (2 ** 2) + Cr * 2 + Dr

    matriz = [
        [variables_a[fila]["valor"] for fila in "abc"]
        + [variables_a["A"]["valor"]],
        [variables_a[fila]["valor"] for fila in "def"]
        + [variables_a["B"]["valor"]],
        [variables_a[fila]["valor"] for fila in "ghi"]
        + [variables_a["C"]["valor"]],
    ]
    sol_sistema, pasos_gj = resolver_gauss_jordan(matriz)
    x, y, z = sol_sistema
    b_s = x + y + z

    # c_s: Regla del trapecio compuesta para el acumulado de l(t) de 0 a 4
    t_vals = (0, 1, 2, 3, 4)
    l_vals = (3.0, 4.4, 6.9, 10.1, 15.0)
    h = t_vals[1] - t_vals[0]
    c_s = (h / 2) * (l_vals[0] + l_vals[-1] + 2 * sum(l_vals[1:-1]))

    # d_s = 10 * Tm + 3 * Tc
    nodos_m = (180, 320, 470, 650)
    valores_tm = (145, 188, 279, 430)
    Tm = interpolacion_lagrange(410, nodos_m, valores_tm)

    matriz_tc = [
        [2**2, 2, 1, 12.4],
        [5**2, 5, 1, 21.7],
        [8**2, 8, 1, 39.5]
    ]
    sol_tc, pasos_gj_tc = resolver_gauss_jordan(matriz_tc)
    P, Q, S_coef = sol_tc
    a_val = variables_a["a"]["valor"]
    Tc = P * (a_val ** 2) + Q * a_val + S_coef

    d_s = 10 * Tm + 3 * Tc

    return {
        "Ar": Ar, "Br": Br, "Cr": Cr, "Dr": Dr,
        "a_s": a_s, "x": x, "y": y, "z": z, "b_s": b_s,
        "c_s": c_s, "Tm": Tm, "Tc": Tc, "d_s": d_s,
        "pasos_gj": pasos_gj, "pasos_gj_tc": pasos_gj_tc,
        "historial_ar": coeficientes["historial_ar"]
    }


def calcular_documento_d():
    """Calcula k (pendiente global por mínimos cuadrados) para el Documento D."""
    r_vals = (0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0)
    y_vals = (1.1, 2.0, 2.4, 3.2, 3.8, 4.1, 5.0, 5.4)
    n = len(r_vals)

    sum_r = sum(r_vals)
    sum_y = sum(y_vals)
    sum_r2 = sum(r * r for r in r_vals)
    sum_ry = sum(r * y for r, y in zip(r_vals, y_vals))

    denom = n * sum_r2 - (sum_r ** 2)
    if abs(denom) < 1e-12:
        raise ValueError("El denominador para el cálculo de k es cero.")

    k = (n * sum_ry - sum_r * sum_y) / denom
    a0 = (sum_y - k * sum_r) / n

    return {
        "r_vals": r_vals, "y_vals": y_vals, "n": n,
        "sum_r": sum_r, "sum_y": sum_y, "sum_r2": sum_r2,
        "sum_ry": sum_ry, "k": k, "a0": a0,
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
        self.title("Sistema de Métodos Numéricos - PIA")
        self.geometry("1100x750")
        self.minsize(850, 620)
        
        # Paleta de colores en Modo Oscuro Sofisticado (No cansa la vista)
        self.colores = {
            "bg_principal": "#0b0f19",       # Fondo general de la app (azul muy oscuro)
            "bg_tarjeta": "#111827",         # Fondo de tarjetas y paneles
            "texto_oscuro": "#f3f4f6",       # Texto principal claro
            "texto_secundario": "#9ca3af",   # Texto secundario gris suave
            "primario": "#3b82f6",           # Azul eléctrico brillante para acentos
            "primario_hover": "#2563eb",     # Azul hover
            "acento": "#06b6d4",             # Cian elegante
            "barra_superior": "#030712",     # Barra inferior/superior oscura
            "borde": "#374151",              # Bordes sutiles
            "linea_tabla": "#1f2937"         # Filas alternadas en tablas
        }
        
        self.configure(bg=self.colores["bg_principal"])
        self.matriculas = {}
        self.resultados = None
        self.entradas = []
        self.mostrar_captura()

    def limpiar_pantalla(self):
        for widget in self.winfo_children():
            widget.destroy()

    def crear_encabezado(self, titulo, subtitulo=""):
        header_frame = tk.Frame(self, bg=self.colores["bg_principal"], padx=35, pady=25)
        header_frame.pack(fill="x")
        
        lbl_titulo = tk.Label(
            header_frame, text=titulo,
            font=("Segoe UI", 22, "bold"), fg=self.colores["texto_oscuro"],
            bg=self.colores["bg_principal"]
        )
        lbl_titulo.pack(anchor="w")
        
        if subtitulo:
            lbl_sub = tk.Label(
                header_frame, text=subtitulo,
                font=("Segoe UI", 11), fg=self.colores["texto_secundario"],
                bg=self.colores["bg_principal"]
            )
            lbl_sub.pack(anchor="w", pady=(2, 0))

    def agregar_navegacion(self):
        barra = tk.Frame(self, bg=self.colores["barra_superior"], height=60)
        barra.pack(side="bottom", fill="x")
        barra.pack_propagate(False)

        botones_nav = [
            ("Resultados", self.mostrar_resultados),
            ("Documento A", self.mostrar_documento_a),
            ("Documento B", self.mostrar_documento_b),
            ("Documento C", self.mostrar_documento_c),
            ("Documento D", self.mostrar_documento_d),
            ("Operaciones detalladas", self.mostrar_operaciones),
        ]

        left_container = tk.Frame(barra, bg=self.colores["barra_superior"])
        left_container.pack(side="left", padx=20)

        for texto, comando in botones_nav:
            btn = tk.Button(
                left_container, text=texto, command=comando,
                bg=self.colores["barra_superior"], fg="#e5e7eb",
                activebackground="#1f2937", activeforeground="white",
                relief="flat", padx=12, pady=8,
                font=("Segoe UI", 10, "bold"), cursor="hand2"
            )
            btn.pack(side="left", padx=4, pady=10)

        right_container = tk.Frame(barra, bg=self.colores["barra_superior"])
        right_container.pack(side="right", padx=20)

        btn_nueva = tk.Button(
            right_container, text="Nueva captura", command=self.mostrar_captura,
            bg=self.colores["acento"], fg="#030712",
            activebackground="#22d3ee", activeforeground="#030712",
            relief="flat", padx=14, pady=7,
            font=("Segoe UI", 10, "bold"), cursor="hand2"
        )
        btn_nueva.pack(pady=12)

    def mostrar_captura(self):
        self.limpiar_pantalla()
        self.configure(bg=self.colores["bg_principal"])

        # Contenedor central principal estilo tarjeta elegante en modo oscuro
        marco_central = tk.Frame(
            self, bg=self.colores["bg_tarjeta"],
            bd=1, relief="solid"
        )
        marco_central.place(relx=0.5, rely=0.5, anchor="center", width=520, height=440)
        
        inner_frame = tk.Frame(marco_central, bg=self.colores["bg_tarjeta"], padx=40, pady=35)
        inner_frame.pack(fill="both", expand=True)

        tk.Label(
            inner_frame, text="Proyecto Integrador de Aplicación",
            font=("Segoe UI", 12, "bold"), fg=self.colores["primario"],
            bg=self.colores["bg_tarjeta"]
        ).pack(anchor="w")

        tk.Label(
            inner_frame, text="Captura de Matrículas",
            font=("Segoe UI", 20, "bold"), fg=self.colores["texto_oscuro"],
            bg=self.colores["bg_tarjeta"]
        ).pack(anchor="w", pady=(5, 8))

        tk.Label(
            inner_frame, text="Ingrese 4 matrículas únicas de 7 dígitos para inicializar los modelos numéricos.",
            font=("Segoe UI", 10), fg=self.colores["texto_secundario"],
            bg=self.colores["bg_tarjeta"], wraplength=420, justify="left"
        ).pack(anchor="w", pady=(0, 25))

        self.entradas = []
        for i in range(1, 5):
            fila_mat = tk.Frame(inner_frame, bg=self.colores["bg_tarjeta"])
            fila_mat.pack(fill="x", pady=6)
            
            lbl = tk.Label(
                fila_mat, text=f"Matrícula M{i}:", width=14, anchor="w",
                font=("Segoe UI", 11, "bold"), fg=self.colores["texto_oscuro"],
                bg=self.colores["bg_tarjeta"]
            )
            lbl.pack(side="left")

            ent = tk.Entry(
                fila_mat, font=("Consolas", 13), width=18,
                justify="center", relief="solid", bd=1,
                bg="#1f2937", fg="#f3f4f6", insertbackground="white"
            )
            ent.pack(side="right", ipady=3)
            self.entradas.append(ent)

        btn_calcular = tk.Button(
            inner_frame, text="Calcular Modelos Numéricos", command=self.procesar,
            bg=self.colores["primario"], fg="white",
            activebackground=self.colores["primario_hover"], activeforeground="white",
            relief="flat", font=("Segoe UI", 11, "bold"), cursor="hand2",
            padx=20, pady=10
        )
        btn_calcular.pack(fill="x", pady=(25, 0))

    def procesar(self):
        valores = [e.get().strip() for e in self.entradas]
        if any(not v.isdigit() or len(v) != 7 for v in valores):
            messagebox.showerror(
                "Error de Validación",
                "Cada matrícula debe contener exactamente 7 dígitos numéricos válidos."
            )
            return
        if len(set(valores)) != 4:
            messagebox.showerror("Error de Validación", "Las cuatro matrículas deben ser distintas.")
            return
        self.matriculas = {
            f"M{i}": [int(d) for d in val]
            for i, val in enumerate(valores, start=1)
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
        nodos_b = (0, P["M2"], P["M3"], P["M4"])
        valores_b = (1, 4, 9, 16)
        B = interpolacion_lagrange(N["M2"], nodos_b, valores_b)
        C, _, _ = resolver_c_newton(P["M1"], P["M2"], S["M3"])
        
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
            "B": {"formula": "Interpolación Lagrange", "valor": B},
            "C": {"formula": "Newton-Raphson u^3 - ... = 0", "valor": C},
        }
        return variables

    def crear_tabla_resultados(self, contenedor, filas_datos, encabezados):
        canvas_tabla = tk.Canvas(contenedor, bg=self.colores["bg_tarjeta"], highlightthickness=0)
        scrollbar = tk.Scrollbar(contenedor, orient="vertical", command=canvas_tabla.yview)
        scrollable_frame = tk.Frame(canvas_tabla, bg=self.colores["bg_tarjeta"])

        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas_tabla.configure(scrollregion=canvas_tabla.bbox("all"))
        )
        canvas_tabla.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas_tabla.configure(yscrollcommand=scrollbar.set)

        # Encabezados
        header_row = tk.Frame(scrollable_frame, bg="#1e3a8a", padx=15, pady=10)
        header_row.pack(fill="x", padx=20, pady=(10, 5))
        for enc in encabezados:
            tk.Label(
                header_row, text=enc, font=("Segoe UI", 11, "bold"),
                fg="#f3f4f6", bg="#1e3a8a", anchor="w", width=25
            ).pack(side="left", padx=10)

        # Filas alternadas oscuras
        for idx, fila in enumerate(filas_datos):
            bg_color = self.colores["bg_tarjeta"] if idx % 2 == 0 else self.colores["linea_tabla"]
            row_frame = tk.Frame(scrollable_frame, bg=bg_color, padx=15, pady=12)
            row_frame.pack(fill="x", padx=20, pady=2)
            
            for val in fila:
                tk.Label(
                    row_frame, text=str(val), font=("Consolas", 11),
                    fg=self.colores["texto_oscuro"], bg=bg_color, anchor="w", width=25
                ).pack(side="left", padx=10)

        canvas_tabla.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

    def mostrar_resultados(self):
        if not self.resultados:
            self.mostrar_captura()
            return
        self.limpiar_pantalla()
        self.configure(bg=self.colores["bg_principal"])
        self.crear_encabezado("Resultados Base", "Parámetros M, N, P y S calculados a partir de las matrículas.")

        contenedor = tk.Frame(self, bg=self.colores["bg_tarjeta"], bd=1, relief="solid")
        contenedor.place(relx=0.05, rely=0.18, relwidth=0.9, relheight=0.68)

        M, N, P, S = self.resultados
        encabezados = ("Matrícula", "M", "N", "P", "S")
        filas = []
        for clave in ("M1", "M2", "M3", "M4"):
            mat_str = "".join(str(d) for d in self.matriculas[clave])
            filas.append((mat_str, str(M[clave]), f"{N[clave]:.7f}", f"{P[clave]:.7f}", str(S[clave])))

        self.crear_tabla_resultados(contenedor, filas, encabezados)
        self.agregar_navegacion()

    def mostrar_documento_a(self):
        if not self.resultados:
            self.mostrar_captura()
            return
        self.limpiar_pantalla()
        self.configure(bg=self.colores["bg_principal"])
        self.crear_encabezado("Documento A: Estado Operativo", "Variables de cálculo y funciones operativas evaluadas.")

        contenedor = tk.Frame(self, bg=self.colores["bg_tarjeta"], bd=1, relief="solid")
        contenedor.place(relx=0.05, rely=0.18, relwidth=0.9, relheight=0.68)

        vars_a = self.calcular_documento_a()
        filas = [(k, f"{v['valor']:.7f}") for k, v in vars_a.items()]
        self.crear_tabla_resultados(contenedor, filas, ("Variable", "Valor numérico"))
        self.agregar_navegacion()

    def mostrar_documento_b(self):
        self.limpiar_pantalla()
        self.configure(bg=self.colores["bg_principal"])
        self.crear_encabezado("Documento B: Respuesta bajo carga", "Coeficientes Ar, Br, Cr y Dr del modelo de respuesta.")

        contenedor = tk.Frame(self, bg=self.colores["bg_tarjeta"], bd=1, relief="solid")
        contenedor.place(relx=0.05, rely=0.18, relwidth=0.9, relheight=0.68)

        coefs_b = calcular_coeficientes_documento_b()
        filas = [
            ("Ar", f"{coefs_b['Ar']:.10f}"),
            ("Br", f"{coefs_b['Br']:.10f}"),
            ("Cr", f"{coefs_b['Cr']:.10f}"),
            ("Dr", f"{coefs_b['Dr']:.10f}"),
        ]
        self.crear_tabla_resultados(contenedor, filas, ("Coeficiente", "Valor calculado"))
        self.agregar_navegacion()

    def mostrar_documento_c(self):
        if not self.resultados:
            self.mostrar_captura()
            return
        self.limpiar_pantalla()
        self.configure(bg=self.colores["bg_principal"])
        self.crear_encabezado("Documento C: Saturación", "Parámetros a_s, b_s, c_s y d_s del sistema de saturación.")

        contenedor = tk.Frame(self, bg=self.colores["bg_tarjeta"], bd=1, relief="solid")
        contenedor.place(relx=0.05, rely=0.18, relwidth=0.9, relheight=0.68)

        vars_a = self.calcular_documento_a()
        val_c = calcular_documento_c(vars_a)
        filas = [
            ("a_s", f"{val_c['a_s']:.7f}"),
            ("b_s", f"{val_c['b_s']:.7f}"),
            ("c_s", f"{val_c['c_s']:.7f}"),
            ("d_s", f"{val_c['d_s']:.7f}"),
        ]
        self.crear_tabla_resultados(contenedor, filas, ("Variable", "Valor final"))
        self.agregar_navegacion()

    def mostrar_documento_d(self):
        self.limpiar_pantalla()
        self.configure(bg=self.colores["bg_principal"])
        self.crear_encabezado("Documento D: Transitorio y EDO", "Parámetro k obtenido mediante ajuste por mínimos cuadrados.")

        contenedor = tk.Frame(self, bg=self.colores["bg_tarjeta"], bd=1, relief="solid")
        contenedor.place(relx=0.05, rely=0.18, relwidth=0.9, relheight=0.68)

        datos_d = calcular_documento_d()
        filas = [("k (Pendiente global)", f"{datos_d['k']:.7f}")]
        self.crear_tabla_resultados(contenedor, filas, ("Parámetro", "Valor"))
        self.agregar_navegacion()

    def mostrar_operaciones(self):
        if not self.resultados:
            self.mostrar_captura()
            return
        self.limpiar_pantalla()
        self.configure(bg=self.colores["bg_principal"])
        self.crear_encabezado("Operaciones e Iteraciones Detalladas", "Desglose completo de matrices, métodos numéricos y fórmulas paso a paso.")

        contenedor = tk.Frame(self, bg=self.colores["bg_tarjeta"], bd=1, relief="solid")
        contenedor.place(relx=0.03, rely=0.16, relwidth=0.94, relheight=0.72)

        texto_area = tk.Text(
            contenedor, bg=self.colores["bg_tarjeta"], fg=self.colores["texto_oscuro"],
            font=("Consolas", 10), relief="flat", padx=20, pady=20, wrap="word",
            insertbackground="white"
        )
        scrollbar = tk.Scrollbar(contenedor, command=texto_area.yview)
        texto_area.configure(yscrollcommand=scrollbar.set)
        texto_area.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

        M, N, P, S = self.resultados
        vars_a = self.calcular_documento_a()
        coefs_b = calcular_coeficientes_documento_b()
        val_c = calcular_documento_c(vars_a)
        datos_d = calcular_documento_d()

        lineas = []
        lineas.append("======================================================================")
        lineas.append(" 1. EXTRACCIÓN Y CÁLCULO DE PARÁMETROS BASE DESDE MATRÍCULAS")
        lineas.append("======================================================================")
        for clave in ("M1", "M2", "M3", "M4"):
            digs = self.matriculas[clave]
            lineas.append(f"• Matrícula {clave} = {''.join(str(d) for d in digs)}")
            lineas.append(f"  - M = {M[clave]}")
            lineas.append(f"  - N = M / 10^7 = {N[clave]:.7f}")
            lineas.append(f"  - P = (1/7) * suma(dígitos) = {P[clave]:.7f}")
            lineas.append(f"  - S = suma(pos * dígito) = {S[clave]}")
            lineas.append("")

        lineas.append("======================================================================")
        lineas.append(" 2. DOCUMENTO A: ESTADO OPERATIVO (CÁLCULO DE VARIABLES Y NEWTON-RAPHSON)")
        lineas.append("======================================================================")
        for k, v in vars_a.items():
            if k == 'C':
                _, _, hist_c = resolver_c_newton(P["M1"], P["M2"], S["M3"])
                lineas.append(f"• Variable C (Raíz por Newton-Raphson de F(u) = u^3 - (P1+P2)u^2 + S3u - 5 = 0):")
                lineas.append(f"  Detalle de iteraciones:")
                for iter_info in hist_c:
                    lineas.append(f"    - Iter {iter_info['iter']}: x_k = {iter_info['x']:.7f}, F(x_k) = {iter_info['fx']:.7f}, F'(x_k) = {iter_info['dfx']:.7f} => x_new = {iter_info['x_nueva']:.7f}")
                lineas.append(f"  -> Resultado final C = {v['valor']:.7f}\n")
            else:
                lineas.append(f"• Variable {k} [{v['formula']}]: {v['valor']:.7f}")
        lineas.append("")

        lineas.append("======================================================================")
        lineas.append(" 3. DOCUMENTO B: RESPUESTA BAJO CARGA (NEWTON-RAPHSON Y MÍNIMOS CUADRADOS)")
        lineas.append("======================================================================")
        lineas.append("• Cálculo de Ar por Newton-Raphson en w^3 - 2.1w^2 + 0.9w - 0.08 = 0:")
        for iter_info in coefs_b['historial_ar']:
            lineas.append(f"  - Iter {iter_info['iter']}: w = {iter_info['x']:.7f}, F(w) = {iter_info['fx']:.7f}, F'(w) = {iter_info['dfx']:.7f}")
        lineas.append(f"  -> Ar final = {coefs_b['Ar']:.10f}\n")
        lineas.append(f"• Cálculo de Br (Integral analítica de 1 a 3): Br = 2*ln(3) + 13/3 = {coefs_b['Br']:.10f}")
        lineas.append(f"• Tendencia lineal global (Cr y Dr) con puntos tabulados:")
        lineas.append(f"  -> Pendiente Cr = {coefs_b['Cr']:.10f}")
        lineas.append(f"  -> Ordenada Dr = {coefs_b['Dr']:.10f}\n")

        lineas.append("======================================================================")
        lineas.append(" 4. DOCUMENTO C: SATURACIÓN (GAUSS-JORDAN PASO A PASO Y TRAPECIO)")
        lineas.append("======================================================================")
        lineas.append(f"• a_s = R(2) = {val_c['a_s']:.7f}")
        lineas.append(f"• b_s = x + y + z (Resolución de Sistema 3x3 por Gauss-Jordan con Pivoteo Parcial):")
        lineas.append(f"  Pasos de reducción de matriz:")
        for paso in val_c['pasos_gj']:
            lineas.append(f"  - [{paso['etapa']}]")
            for fila in paso['matriz']:
                lineas.append(f"    { [round(val, 4) for val in fila] }")
        lineas.append(f"  -> Solución del sistema: x = {val_c['x']:.5f}, y = {val_c['y']:.5f}, z = {val_c['z']:.5f}")
        lineas.append(f"  -> b_s = x + y + z = {val_c['b_s']:.7f}\n")
        lineas.append(f"• c_s (Acumulado por Regla del Trapecio Compuesta de l(t) de 0 a 4): {val_c['c_s']:.7f}")
        lineas.append(f"• d_s = 10 * Tm + 3 * Tc:")
        lineas.append(f"  - Tm (Interpolación Lagrange en m=410) = {val_c['Tm']:.7f}")
        lineas.append(f"  - Tc (Evaluación cuadrática con lambda = a) = {val_c['Tc']:.7f}")
        lineas.append(f"  -> d_s final = {val_c['d_s']:.7f}\n")

        lineas.append("======================================================================")
        lineas.append(" 5. DOCUMENTO D: TRANSITORIO Y PENDIENTE GLOBAL")
        lineas.append("======================================================================")
        lineas.append("• Cálculo de k (Pendiente global por Mínimos Cuadrados):")
        lineas.append(f"  - Número de puntos n = {datos_d['n']}")
        lineas.append(f"  - Sumatorias: ∑r = {datos_d['sum_r']:.2f}, ∑y = {datos_d['sum_y']:.2f}, ∑r² = {datos_d['sum_r2']:.2f}, ∑ry = {datos_d['sum_ry']:.2f}")
        lineas.append(f"  - Sustitución: k = [n*(∑ry) - (∑r)*(∑y)] / [n*(∑r²) - (∑r)²]")
        lineas.append(f"  -> k final = {datos_d['k']:.7f} (Ordenada a0 = {datos_d['a0']:.7f})")

        texto_area.insert("1.0", "\n".join(lineas))
        texto_area.configure(state="disabled")
        self.agregar_navegacion()


if __name__ == "__main__":
    Aplicacion().mainloop()