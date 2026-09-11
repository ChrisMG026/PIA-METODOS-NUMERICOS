import tkinter as tk
from tkinter import messagebox


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
        self.geometry("900x620")
        self.minsize(720, 500)
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

    def agregar_navegacion(self, canvas):
        barra = tk.Frame(self, bg="#172033", height=52)
        barra.place(relx=0, rely=1, relwidth=1, anchor="sw")
        tk.Button(
            barra, text="Resultados", command=self.mostrar_resultados,
            bg="#172033", fg="white", relief="flat", padx=18,
            font=("Segoe UI", 10, "bold")
        ).pack(side="left", padx=18, pady=10)
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