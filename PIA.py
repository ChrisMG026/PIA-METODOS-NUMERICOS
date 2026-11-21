import tkinter as tk
import tkinter.font as tkfont
from tkinter import ttk
import math
from tkinter import filedialog, messagebox
from html import escape
from pathlib import Path

#código de funciones matemáticas y de cálculo de documentos

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
    if len(nodos) != len(valores) or not nodos:
        raise ValueError("La interpolación requiere nodos y valores del mismo tamaño.")
    if len(set(nodos)) != len(nodos):
        raise ValueError("Los nodos de interpolación deben ser distintos.")
    resultado = 0.0
    for indice, nodo in enumerate(nodos):
        termino = valores[indice]
        for otro_indice, otro_nodo in enumerate(nodos):
            if indice != otro_indice:
                termino *= (x - otro_nodo) / (nodo - otro_nodo)
        resultado += termino
    return resultado


def interpolacion_newton_adelante(x, nodos, valores):
    """Interpola con diferencias hacia adelante sobre una malla equidistante."""
    if len(nodos) != len(valores) or not nodos:
        raise ValueError("La interpolación requiere nodos y valores del mismo tamaño.")
    if len(nodos) == 1:
        return valores[0]
    paso = nodos[1] - nodos[0]
    if paso == 0 or any(
        not math.isclose(nodos[indice + 1] - nodos[indice], paso, rel_tol=1e-10, abs_tol=1e-12)
        for indice in range(1, len(nodos) - 1)
    ):
        raise ValueError("Newton hacia adelante requiere nodos equidistantes.")

    diferencias = [list(valores)]
    while len(diferencias[-1]) > 1:
        anterior = diferencias[-1]
        diferencias.append([
            anterior[indice + 1] - anterior[indice]
            for indice in range(len(anterior) - 1)
        ])

    parametro = (x - nodos[0]) / paso
    producto = 1.0
    resultado = diferencias[0][0]
    for orden in range(1, len(diferencias)):
        producto *= (parametro - orden + 1) / orden
        resultado += producto * diferencias[orden][0]
    return resultado


def resolver_c_newton(
    P1, P2, S3, tolerancia=1e-10, tolerancia_residuo=1e-10,
    max_iteraciones=100,
):
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
        if not math.isfinite(nueva_aproximacion):
            raise OverflowError("Newton-Raphson produjo una aproximación no finita.")
        residuo_nuevo = (
            nueva_aproximacion ** 3
            - suma_p * nueva_aproximacion ** 2
            + S3 * nueva_aproximacion
            - 5
        )
        historial.append({
            "iter": iteracion,
            "x": aproximacion,
            "fx": funcion,
            "dfx": derivada,
            "x_nueva": nueva_aproximacion,
            "residuo_nuevo": residuo_nuevo,
        })
        if (
            abs(nueva_aproximacion - aproximacion) < tolerancia
            and abs(residuo_nuevo) < tolerancia_residuo
        ):
            return nueva_aproximacion, iteracion, historial
        aproximacion = nueva_aproximacion
    raise ValueError("Newton-Raphson no convergio en el limite de iteraciones.")


def calcular_coeficientes_documento_b():
    """Calcula Ar, Br, Cr y Dr para el modelo de respuesta R(n)."""
    raices_ar = raices_reales_cubica(1, -2.1, 0.9, -0.08, tolerancia=1e-12)
    raices_positivas_ar = [raiz for raiz in raices_ar if raiz > 0]
    if not raices_positivas_ar:
        raise ValueError("El polinomio de Ar no tiene raíces positivas.")
    Ar = min(raices_positivas_ar)

    # ---------- Br: integral definida de 1 a 3 de (2/x + 0.5x^2) dx ----------
    # Antiderivada: 2*ln|x| + x^3/6  →  evaluada: [2ln3 + 27/6] - [0 + 1/6] = 2ln3 + 13/3
    Br = 2 * ln(3) + 13 / 3

    # ---------- Cr: tendencia lineal global sobre la tabla de mediciones ----------
    nodos_cr = (1, 2, 4, 6)
    valores_cr = (4.2, 6.1, 9.8, 13.9)
    media_n = sum(nodos_cr) / len(nodos_cr)
    media_h = sum(valores_cr) / len(valores_cr)
    Cr = sum(
        (nodo - media_n) * (valor - media_h)
        for nodo, valor in zip(nodos_cr, valores_cr)
    ) / sum((nodo - media_n) ** 2 for nodo in nodos_cr)

    # ---------- Dr: evaluar en x = 2.5 la tabla propia de Dr (interpolación) ----------
    nodos_dr = (1, 2, 3, 4)
    valores_dr = (2.8, 4.6, 7.9, 12.5)
    Dr = interpolacion_newton_adelante(2.5, nodos_dr, valores_dr)

    return {
        "Ar": Ar,
        "Br": Br,
        "Cr": Cr,
        "Dr": Dr,
        "raices_ar": raices_ar,
        "historial_ar": [],
        "nodos_cr": nodos_cr,
        "valores_cr": valores_cr,
        "nodos_dr": nodos_dr,
        "valores_dr": valores_dr,
        "media_n": media_n,
        "media_h": media_h,
    }

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


def calcular_documento_c(variables_a, coeficientes_b=None):
    """Calcula a_s, b_s, c_s y d_s para el Documento C."""
    coeficientes = coeficientes_b or calcular_coeficientes_documento_b()
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

    # c_s: Simpson 1/3 compuesto para el acumulado de l(t) de 0 a 4
    t_vals = (0, 1, 2, 3, 4)
    l_vals = (3.0, 4.4, 6.9, 10.1, 15.0)
    h = t_vals[1] - t_vals[0]
    desglose_c_s = desglose_simpson_1_3_tabla(l_vals, h)
    c_s = desglose_c_s["resultado"]

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
        "desglose_c_s": {**desglose_c_s, "nodos": t_vals, "valores": l_vals},
        "pasos_gj": pasos_gj, "pasos_gj_tc": pasos_gj_tc,
        "historial_ar": coeficientes["historial_ar"]
    }


def calcular_documento_d():
    """Calcula k (pendiente global por mínimos cuadrados) para el Documento D."""
    r_vals = (1, 3, 5, 7)
    y_vals = (2.0, 5.3, 8.9, 12.2)
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
        "r_vals": r_vals, "y_vals": y_vals, "n_puntos": n,
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

# ==========================================
# NUEVOS MÉTODOS NUMÉRICOS (CLASE)
# ==========================================
def metodo_biseccion_historial(func, a, b, tolerancia=1e-7, max_iter=100):
    if tolerancia <= 0 or max_iter <= 0:
        raise ValueError("La tolerancia y el máximo de iteraciones deben ser positivos.")
    if a >= b:
        raise ValueError("El intervalo de bisección debe cumplir a < b.")
    fa, fb = func(a), func(b)
    if fa == 0:
        return a, [{"iter": 0, "a": a, "b": b, "x": a, "fx": fa}]
    if fb == 0:
        return b, [{"iter": 0, "a": a, "b": b, "x": b, "fx": fb}]
    if fa * fb > 0:
        raise ValueError("El intervalo no contiene una raíz.")
    historial = []
    for iteracion in range(1, max_iter + 1):
        c = (a + b) / 2.0
        fc = func(c)
        historial.append({
            "iter": iteracion, "a": a, "b": b, "x": c,
            "fx": fc, "ancho": b-a,
        })
        if fc == 0 or abs(fc) < tolerancia or (b - a) / 2.0 < tolerancia:
            return c, historial
        if fa * fc < 0:
            b = c
            fb = fc
        else:
            a = c
            fa = fc
    raise ValueError("Bisección no convergió en el límite de iteraciones.")

def metodo_biseccion(func, a, b, tolerancia=1e-7, max_iter=100):
    raiz, _ = metodo_biseccion_historial(func, a, b, tolerancia, max_iter)
    return raiz

def simpson_1_3(func, a, b, n):
    if n <= 0 or n % 2 != 0:
        raise ValueError("n debe ser un entero positivo y par para Simpson 1/3.")
    h = (b - a) / n
    suma = func(a) + func(b)
    for i in range(1, n):
        x = a + i * h
        suma += (2 if i % 2 == 0 else 4) * func(x)
    return (h / 3) * suma

def desglose_simpson_1_3_tabla(valores, h):
    """Desglosa Simpson 1/3 compuesto aplicado a una tabla equidistante.

    Devuelve las sumatorias ponderadas y el resultado para poder mostrar todo
    el proceso en la interfaz y en el informe.
    """
    if h <= 0 or len(valores) < 3 or len(valores) % 2 == 0:
        raise ValueError(
            "Simpson 1/3 compuesto requiere al menos 3 valores, una cantidad "
            "impar de valores (par de intervalos) y paso positivo."
        )
    extremos = valores[0] + valores[-1]
    impares = sum(valores[1:-1:2])
    pares = sum(valores[2:-1:2])
    suma_ponderada = extremos + 4 * impares + 2 * pares
    return {
        "h": h,
        "num_intervalos": len(valores) - 1,
        "suma_extremos": extremos,
        "suma_impares": impares,
        "suma_pares": pares,
        "suma_ponderada": suma_ponderada,
        "resultado": (h / 3) * suma_ponderada,
    }


def simpson_1_3_tabla(valores, h):
    """Aplica Simpson 1/3 compuesto a valores tabulados con paso constante."""
    return desglose_simpson_1_3_tabla(valores, h)["resultado"]


def simpson_3_8(valores, h):
    if len(valores) != 4 or h <= 0:
        raise ValueError("Simpson 3/8 requiere cuatro valores y un paso positivo.")
    return (3 * h / 8) * (valores[0] + 3*valores[1] + 3*valores[2] + valores[3])

def euler_adelante_historial(edo, t0, y0, t_final, h):
    if h <= 0 or t_final < t0:
        raise ValueError("El paso debe ser positivo y t_final no puede preceder a t0.")
    t, y = t0, y0
    historial = [(t, y)]
    while t < t_final:
        paso = min(h, t_final - t)
        y = y + paso * edo(t, y)
        if not math.isfinite(y):
            raise OverflowError("Euler produjo un valor no finito; reduzca el paso o revise el modelo.")
        t = min(t_final, t + paso)
        historial.append((t, y))
    return historial


class IntegracionNoFinita(ArithmeticError):
    def __init__(self, tiempo, historial, motivo):
        super().__init__(motivo)
        self.tiempo = tiempo
        self.historial = historial


def runge_kutta_4_historial(edo, t0, y0, t_final, h):
    """Integra una EDO con RK4 clásico y devuelve los pasos aceptados."""
    if h <= 0 or t_final < t0:
        raise ValueError("El paso debe ser positivo y t_final no puede preceder a t0.")
    if not all(math.isfinite(valor) for valor in (t0, y0, t_final, h)):
        raise ValueError("Los datos de RK4 deben ser finitos.")

    historial = [(t0, y0)]
    tiempo = t0

    def evaluar(t, valor):
        try:
            pendiente = edo(t, valor)
        except (OverflowError, ZeroDivisionError) as error:
            raise IntegracionNoFinita(
                t, historial, f"El cálculo de la pendiente falló: {error}"
            ) from error
        if not math.isfinite(pendiente):
            raise IntegracionNoFinita(
                t, historial, "La pendiente de la EDO dejó de ser finita."
            )
        return pendiente

    while tiempo < t_final:
        paso = min(h, t_final - tiempo)
        if tiempo + paso <= tiempo:
            raise IntegracionNoFinita(
                tiempo, historial, "El tamaño del paso ya no permite avanzar en el tiempo."
            )
        try:
            k1 = evaluar(tiempo, historial[-1][1])
            k2 = evaluar(tiempo + paso / 2, historial[-1][1] + paso * k1 / 2)
            k3 = evaluar(tiempo + paso / 2, historial[-1][1] + paso * k2 / 2)
            k4 = evaluar(tiempo + paso, historial[-1][1] + paso * k3)
            siguiente = historial[-1][1] + paso * (k1 + 2*k2 + 2*k3 + k4) / 6
        except OverflowError as error:
            raise IntegracionNoFinita(
                tiempo, historial, "RK4 produjo un desbordamiento numérico."
            ) from error
        if not math.isfinite(siguiente):
            raise IntegracionNoFinita(
                tiempo, historial, "RK4 produjo un estado no finito."
            )
        tiempo = min(t_final, tiempo + paso)
        historial.append((tiempo, siguiente))
    return historial


def euler_adelante(edo, t0, y0, t_final, h):
    return euler_adelante_historial(edo, t0, y0, t_final, h)[-1][1]

def interpolar_historial(historial, t):
    if not historial or t < historial[0][0] or t > historial[-1][0]:
        raise ValueError("El tiempo solicitado está fuera del historial.")
    for (t0, y0), (t1, y1) in zip(historial, historial[1:]):
        if t0 <= t <= t1:
            if t1 == t0:
                return y0
            return y0 + (y1 - y0) * (t - t0) / (t1 - t0)
    return historial[-1][1]

def raices_reales_cubica(a, b, c, d, tolerancia=1e-10):
    if abs(a) < 1e-14:
        raise ValueError("El coeficiente cúbico no puede ser cero.")

    def funcion(x):
        return ((a*x + b)*x + c)*x + d

    limite = 1 + max(abs(b/a), abs(c/a), abs(d/a))
    discriminante = 4*b*b - 12*a*c
    puntos_criticos = []
    if discriminante >= 0:
        raiz_discriminante = discriminante ** 0.5
        puntos_criticos = [
            (-2*b - raiz_discriminante) / (6*a),
            (-2*b + raiz_discriminante) / (6*a),
        ]
    puntos = [-limite] + sorted(
        punto for punto in puntos_criticos if -limite < punto < limite
    ) + [limite]
    raices = []

    for punto in puntos:
        if abs(funcion(punto)) <= tolerancia * max(1.0, abs(a*punto**3), abs(d)):
            raices.append(punto)

    for inicio, fin in zip(puntos, puntos[1:]):
        if funcion(inicio) * funcion(fin) < 0:
            raices.append(metodo_biseccion(funcion, inicio, fin, tolerancia))

    raices.sort()
    unicas = []
    for raiz in raices:
        if not unicas or abs(raiz - unicas[-1]) > tolerancia * 10:
            unicas.append(raiz)
    return unicas


def raices_reales_cuadratica(a, b, c):
    if abs(a) < 1e-14:
        raise ValueError("El coeficiente cuadrático no puede ser cero.")
    discriminante = b*b - 4*a*c
    if discriminante < 0:
        return []
    raiz_discriminante = math.sqrt(discriminante)
    if raiz_discriminante == 0:
        return [-b / (2*a)]
    q = -0.5 * (b + math.copysign(raiz_discriminante, b))
    return sorted((q / a, c / q))

# ==========================================
# NUEVOS DOCUMENTOS (D, E, F, G)
# ==========================================
def calcular_documento_d_completo():
    datos_regresion = calcular_documento_d()
    m_val = simpson_1_3(lambda u: 3*u**2 + 1.5*u + 2, 0, 2, 2)
    n_val, historial_n = metodo_biseccion_historial(
        lambda x: x**3 - 7*x + 6, 1 + 1e-8, 3, tolerancia=1e-10
    )
    historial_p = runge_kutta_4_historial(
        lambda t, y: -0.6*y + 4*t, 0, 1, 1.5, 0.1
    )
    p_val = historial_p[-1][1]
    q_val = interpolacion_newton_adelante(
        4.5, (2, 4, 6, 8), (10.5, 18.7, 31.1, 47.3)
    )

    datos_balance = {
        **datos_regresion, "m": m_val, "n": n_val, "p": p_val, "q": q_val,
        "historial_n": historial_n, "iteraciones_n": len(historial_n),
        "historial_p": historial_p,
    }
    datos_balance["L3"] = (
        datos_balance["k"] * 3**4
        + datos_balance["m"] * 3**3
        + datos_balance["n"] * 3**2
        + datos_balance["p"] * 3
        + datos_balance["q"]
    )
    return datos_balance

def calcular_documento_e_completo(val_c, datos_d=None, h=0.1):
    if datos_d is None:
        datos_d = calcular_documento_d_completo()
    Ar, Tm, x = val_c["Ar"], val_c["Tm"], val_c["x"]

    raices_sigma = raices_reales_cuadratica(1, -5.2, 4.8)
    raices_positivas = [raiz for raiz in raices_sigma if raiz > 0]
    if not raices_positivas:
        raise ValueError("La ecuación de r_sigma no tiene raíces positivas.")
    r_sigma = min(raices_positivas)
    sigma = r_sigma / 10.0
    tau = Tm - Ar
    kappa = datos_d["L3"]
    omega = simpson_3_8((1.5, 2.6, 4.8, 8.3), 1.0)

    R_1_8 = Ar*(1.8**3) + val_c["Br"]*(1.8**2) + val_c["Cr"]*1.8 + val_c["Dr"]
    Y0 = R_1_8 + x

    def edo(t, Y):
        if t < 0:
            raise ValueError("El tiempo del modelo E debe ser no negativo.")
        return -sigma*(Y - tau) + kappa/omega

    historial = runge_kutta_4_historial(edo, 0, Y0, 3.0, h)

    return {
        "r_sigma": r_sigma, "raices_sigma": raices_sigma,
        "sigma": sigma, "tau": tau,
        "kappa": kappa, "omega": omega, "Y0": Y0,
        "Y_final": historial[-1][1], "historial": historial,
        "h": h,
        "coeficientes_saturacion": {
            "a_s": val_c["a_s"], "b_s": val_c["b_s"],
            "c_s": val_c["c_s"], "d_s": val_c["d_s"],
        },
    }

def calcular_documento_f_completo(variables_a, val_e, h=0.1):
    u_vals, k_vals = (0, 1, 2, 3), (0.8, 1.1, 1.5, 1.9)
    u_media = sum(u_vals) / len(u_vals)
    k_media = sum(k_vals) / len(k_vals)
    chi = sum((u-u_media)*(k-k_media) for u, k in zip(u_vals, k_vals)) / sum(
        (u-u_media)**2 for u in u_vals
    )
    raices_psi = raices_reales_cuadratica(1, -11, 24)
    if not raices_psi:
        raise ValueError("La ecuación de psi no tiene raíces reales.")
    psi = max(raices_psi)
    xi = variables_a["a"]["valor"] + variables_a["d"]["valor"] + variables_a["g"]["valor"]
    if abs(xi) < 1e-14:
        raise ValueError("xi = a + d + g no puede ser cero.")

    theta0 = 25.0
    capacidad_maxima = chi * psi**2 / 4
    fuentes = [
        interpolar_historial(val_e["historial"], t) / xi
        for t, _ in val_e["historial"]
    ]
    fuente_minima, fuente_maxima = min(fuentes), max(fuentes)
    edo = lambda t, T: (
        chi*T*(psi - T) - interpolar_historial(val_e["historial"], t)/xi
    )
    causa = None
    tiempo_fallo = None
    try:
        historial = runge_kutta_4_historial(edo, 0, theta0, 3.0, h)
        estado = "completo"
    except IntegracionNoFinita as error:
        historial = error.historial
        tiempo_fallo = error.tiempo
        detalle_error = str(error).rstrip(". ")
        if fuente_minima > capacidad_maxima:
            causa = (
                f"{detalle_error}. El aporte Y/xi se mantiene por encima de la "
                f"capacidad máxima del término logístico "
                f"({fuente_minima:.8g} > {capacidad_maxima:.8g}) en [0, 3]. "
                f"Por la ecuación original, theta' <= "
                f"{capacidad_maxima-fuente_minima:.8g} < 0; al caer theta "
                "por debajo de cero, el término "
                "logístico también se vuelve negativo y la magnitud crece "
                "hasta que RK4 deja de obtener valores finitos."
            )
        elif fuente_maxima > capacidad_maxima:
            causa = (
                f"{detalle_error}. El aporte Y/xi alcanza {fuente_maxima:.8g}, por "
                f"encima de la capacidad máxima del término logístico "
                f"({capacidad_maxima:.8g}) en parte del intervalo; esto "
                "contribuye a que la derivada neta sea negativa."
            )
        else:
            causa = str(error)
        estado = "divergente"

    return {
        "estado": estado, "causa": causa, "tiempo_fallo": tiempo_fallo,
        "chi": chi, "psi": psi, "raices_psi": raices_psi,
        "xi": xi, "theta0": theta0,
        "theta_final": historial[-1][1] if estado == "completo" else None,
        "historial": historial, "h": h,
        "fuente_minima": fuente_minima,
        "fuente_maxima": fuente_maxima,
        "capacidad_maxima": capacidad_maxima,
        "justificacion_theta0": (
            "Se propone theta(0)=25 °C como temperatura ambiente/de reposo "
            "típica antes de iniciar la carga de trabajo."
        ),
    }

def calcular_documento_g_completo(val_e, val_f):
    if val_f["estado"] != "completo":
        return {
            "estado": "no disponible",
            "causa": (
                "El modelo F original no produjo una trayectoria finita en "
                f"[0, 3]; la integración se detuvo cerca de t={val_f['tiempo_fallo']:.6g}. "
                "No es válido calcular theta(1) ni el índice G."
            ),
        }
    coeficientes = val_e["coeficientes_saturacion"]
    raices = raices_reales_cubica(
        coeficientes["a_s"], -coeficientes["b_s"],
        coeficientes["c_s"], -coeficientes["d_s"]
    )
    raices_positivas = [raiz for raiz in raices if raiz > 0]
    if not raices_positivas:
        raise ValueError("El modelo de saturación no tiene raíces reales positivas.")
    U_star = min(raices_positivas)
    intervalo = 2.0
    n = round(intervalo / val_e["h"])
    if n % 2:
        n += 1
    integral_Y = simpson_1_3(
        lambda t: interpolar_historial(val_e["historial"], t),
        0, intervalo, n
    )
    theta_1 = interpolar_historial(val_f["historial"], 1.0)
    S = (1/2)*integral_Y + (1/3)*theta_1 + ln(1 + U_star)
    return {
        "estado": "completo",
        "integral_Y": integral_Y, "theta_1": theta_1,
        "U_star": U_star, "raices_saturacion": raices,
        "S": S,
    }

def calcular_sensibilidad(
    resultados, variables_a, val_c, datos_b, datos_d,
    pasos=(0.1, 0.05, 0.025),
):
    sensibilidad_pasos = []
    for h in pasos:
        val_e = calcular_documento_e_completo(val_c, datos_d, h)
        val_f = calcular_documento_f_completo(variables_a, val_e, h)
        val_g = calcular_documento_g_completo(val_e, val_f)
        sensibilidad_pasos.append({
            "h": h, "Y_final": val_e["Y_final"],
            "estado_F": val_f["estado"],
            "tiempo_fallo_F": val_f["tiempo_fallo"],
            "theta_final": val_f["theta_final"],
            "S": val_g.get("S"),
            "estado_G": val_g["estado"],
        })
    for indice, resultado in enumerate(sensibilidad_pasos):
        siguiente = sensibilidad_pasos[indice + 1] if indice + 1 < len(sensibilidad_pasos) else None
        resultado["diferencia_S"] = (
            abs(resultado["S"] - siguiente["S"])
            if siguiente and resultado["S"] is not None and siguiente["S"] is not None
            else None
        )
        resultado["diferencia_Y"] = (
            abs(resultado["Y_final"] - siguiente["Y_final"])
            if siguiente else None
        )

    P, S_base = resultados[2], resultados[3]
    sensibilidad_tolerancias = []
    for tolerancia in (1e-6, 1e-8, 1e-10):
        variables_tolerancia = calcular_variables_documento_a(
            resultados, tolerancia_newton=tolerancia
        )
        val_c_tolerancia = calcular_documento_c(variables_tolerancia, datos_b)
        val_e_tolerancia = calcular_documento_e_completo(
            val_c_tolerancia, datos_d
        )
        val_f_tolerancia = calcular_documento_f_completo(
            variables_tolerancia, val_e_tolerancia
        )
        val_g_tolerancia = calcular_documento_g_completo(
            val_e_tolerancia, val_f_tolerancia
        )
        raiz = variables_tolerancia["C"]["valor"]
        historial = variables_tolerancia["C"]["historial"]
        residuo = (
            raiz**3 - (P["M1"] + P["M2"])*raiz**2
            + S_base["M3"]*raiz - 5
        )
        sensibilidad_tolerancias.append({
            "tolerancia": tolerancia,
            "C": raiz,
            "iteraciones": len(historial),
            "residuo": residuo,
            "estado_F": val_f_tolerancia["estado"],
            "S": val_g_tolerancia.get("S"),
            "diferencia_S": (
                abs(val_g_tolerancia["S"] - val_g["S"])
                if val_g and val_g.get("estado") == "completo"
                and val_g_tolerancia["estado"] == "completo"
                else None
            ),
        })
    return {"pasos": sensibilidad_pasos, "tolerancias": sensibilidad_tolerancias}

def calcular_variables_documento_a(resultados, tolerancia_newton=1e-10):
    _, N, P, S = resultados
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
    C, iteraciones_c, historial_c = resolver_c_newton(
        P["M1"], P["M2"], S["M3"], tolerancia=tolerancia_newton
    )

    return {
        "a": {"formula": "P1 + N2", "valor": a},
        "b": {"formula": "sqrt(S2)", "valor": b},
        "c": {"formula": "P3^2 / P4", "valor": c},
        "d": {"formula": "ln(1 + S1)", "valor": d},
        "e": {"formula": "P2^2 + N3", "valor": e},
        "f": {"formula": "S4 / 10", "valor": f},
        "g": {"formula": "exp(N1)", "valor": g},
        "h": {"formula": "S3 / (P1 + 1)", "valor": h},
        "i": {"formula": "sqrt(P4 + S2)", "valor": i},
        "A": {"formula": "Integral de 0 a P1 de (2t^2 + S2t + 1)dt", "valor": A},
        "B": {"formula": "Interpolación Lagrange en N2", "valor": B},
        "C": {
            "formula": "Newton-Raphson", "valor": C,
            "iteraciones": iteraciones_c, "historial": historial_c,
        },
    }

def calcular_modelo(matriculas):
    resultados = calcular_resultados(matriculas)
    variables_a = calcular_variables_documento_a(resultados)
    datos_b = calcular_coeficientes_documento_b()
    val_c = calcular_documento_c(variables_a, datos_b)
    datos_d = calcular_documento_d_completo()
    val_e = calcular_documento_e_completo(val_c, datos_d)
    val_f = calcular_documento_f_completo(variables_a, val_e)
    val_g = calcular_documento_g_completo(val_e, val_f)
    sensibilidad = calcular_sensibilidad(
        resultados, variables_a, val_c, datos_b, datos_d
    )
    return {
        "resultados": resultados,
        "variables_a": variables_a,
        "documento_b": datos_b,
        "documento_c": val_c,
        "documento_d": datos_d,
        "documento_e": val_e,
        "documento_f": val_f,
        "documento_g": val_g,
        "sensibilidad": sensibilidad,
    }

SIGNIFICADOS_BASE = {
    "M": "Los siete dígitos de la matrícula escritos como un solo número.",
    "N": "La matrícula expresada en una escala pequeña para usarla en los modelos.",
    "P": "El promedio de los siete dígitos de la matrícula.",
    "S": "La suma de los dígitos dando más peso a los que están al final.",
}

SIGNIFICADOS_A = {
    "a": "Combina el promedio de la primera matrícula con la escala de la segunda.",
    "b": "Raíz cuadrada del puntaje ponderado de la segunda matrícula.",
    "c": "Relación entre los promedios de las matrículas tercera y cuarta.",
    "d": "Logaritmo del puntaje ponderado de la primera matrícula más uno.",
    "e": "Promedio de la segunda matrícula al cuadrado más la escala de la tercera.",
    "f": "Una décima parte del puntaje ponderado de la cuarta matrícula.",
    "g": "Exponencial de la escala correspondiente a la primera matrícula.",
    "h": "Relación entre el puntaje de la tercera y el promedio de la primera.",
    "i": "Raíz de la suma del promedio de la cuarta y el puntaje de la segunda.",
    "A": "Resultado de integrar la expresión recuperada para construir el sistema.",
    "B": "Estimación obtenida de la tabla de puntos del expediente.",
    "C": "Solución de la ecuación cúbica del estado operativo.",
}

SIGNIFICADOS_B = {
    "Ar": "Tamaño del término cúbico del modelo de respuesta bajo carga.",
    "Br": "Aporte calculado mediante el área bajo una curva entre 1 y 3.",
    "Cr": "Pendiente que resume la tendencia general de las mediciones.",
    "Dr": "Valor estimado a mitad de la tabla de mediciones.",
}

SIGNIFICADOS_C = {
    "a_s": "Respuesta del modelo bajo carga cuando la entrada vale 2.",
    "b_s": "Suma de las tres incógnitas encontradas al resolver el sistema.",
    "c_s": "Acumulado de las mediciones entre 0 y 4 con Simpson 1/3 compuesto.",
    "d_s": "Combinación de las estimaciones Tm y Tc calculadas con los datos del expediente.",
}

SIGNIFICADOS_D = {
    "k": "Cambio promedio de las mediciones por cada unidad de r.",
    "a0": "Valor de y que predice la recta cuando r vale cero.",
    "m": "Área estimada bajo la curva del modelo entre 0 y 2.",
    "n": "Valor entre 1 y 3 que hace que la ecuación se acerque a cero.",
    "p": "Valor estimado de y cuando el tiempo llega a 1.5.",
    "q": "Valor estimado por la tabla cuando r vale 4.5.",
    "L(3)": "Valor del polinomio de balance del Documento D cuando t vale 3.",
}

SIGNIFICADOS_E = {
    "r_sigma": "Raíz elegida para definir el comportamiento de la respuesta.",
    "sigma": "Factor que determina qué tan rápido cambia la respuesta.",
    "tau": "Nivel objetivo calculado con los datos recuperados.",
    "kappa": "Valor del balance acumulado L(3), calculado con los coeficientes del Documento D.",
    "omega": "Acumulado de la tabla v(t) entre 0 y 3, aproximado con Simpson 3/8.",
    "Y(0)": "Respuesta calculada al inicio del intervalo.",
    "Y(3)": "Respuesta estimada después de tres unidades de tiempo.",
}

SIGNIFICADOS_F = {
    "chi": "Qué tan rápido cambia la temperatura según las mediciones disponibles.",
    "psi": "Valor elegido para definir el límite de temperatura del modelo.",
    "xi": "Combinación de tres valores operativos calculados en el Documento A.",
    "Y/xi": "Intensidad de la respuesta que aparece en la ecuación térmica original.",
    "capacidad": "Límite máximo del término logístico de la ecuación térmica.",
    "theta(0)": "Temperatura inicial propuesta: 25 °C antes de iniciar la carga.",
    "theta(3)": "Temperatura estimada al final de las tres unidades de tiempo.",
    "fuente_minima": "Menor valor observado de Y/xi durante el intervalo.",
    "estado": "Indica si el modelo térmico pudo integrarse hasta t=3.",
}

SIGNIFICADOS_G = {
    "Integral de Y": "Aporte acumulado de la respuesta del sistema entre 0 y 2.",
    "theta(1)": "Temperatura estimada cuando ha transcurrido una unidad de tiempo.",
    "U*": "Menor solución positiva de la ecuación de saturación; representa el nivel operativo elegido.",
    "S": "Índice general que combina respuesta, temperatura y saturación.",
}

def crear_filas_descritas(datos, formato=".7f"):
    return [
        (
            nombre,
            valor if isinstance(valor, str) else (
                "No disponible" if valor is None else format(valor, formato)
            ),
            significado,
        )
        for nombre, valor, significado in datos
    ]

def limitar_trayectoria_para_grafica(valores, estado):
    valores_graficados = list(valores)
    if estado != "divergente" or len(valores) < 2:
        return valores_graficados, False

    indice_salto = None
    for indice in range(1, len(valores)):
        previos = valores[:indice]
        escala_previa = max(
            max(abs(valor) for valor in previos),
            max(previos) - min(previos),
            1.0,
        )
        if abs(valores[indice] - valores[indice - 1]) > 20 * escala_previa:
            indice_salto = indice
            break
    if indice_salto is None:
        return valores_graficados, False

    referencia = valores[:indice_salto]
    minimo = min(referencia)
    maximo = max(referencia)
    amplitud = max(maximo - minimo, abs(minimo) * 0.1, abs(maximo) * 0.1, 1e-9)
    limite_inferior = minimo - amplitud * 0.15
    limite_superior = maximo + amplitud * 0.15
    for indice in range(indice_salto, len(valores)):
        valor = valores[indice]
        valor_limitado = min(limite_superior, max(limite_inferior, valor))
        if valor_limitado != valor:
            valores_graficados[indice] = valor_limitado
    return valores_graficados, valores_graficados != valores


def crear_svg_trayectoria(
    historial, titulo, color, eje_y="Valor", estado=None, tiempo_fallo=None
):
    ancho, alto = 720, 360
    izquierda, derecha, arriba, abajo = 76, 24, 48, 62
    x0, x1 = izquierda, ancho - derecha
    y0, y1 = arriba, alto - abajo
    tiempos = [punto[0] for punto in historial]
    valores = [punto[1] for punto in historial]
    t_min, t_max = min(tiempos), max(tiempos)
    valores_graficados, recorte_escala = limitar_trayectoria_para_grafica(
        valores, estado
    )
    y_min, y_max = min(valores_graficados), max(valores_graficados)
    if y_min == y_max:
        y_max = y_min + 1

    def x_pixel(tiempo):
        return x0 + (tiempo - t_min) * (x1 - x0) / max(t_max - t_min, 1e-12)

    def y_pixel(valor):
        return y1 - (valor - y_min) * (y1 - y0) / (y_max - y_min)

    puntos = [
        (x_pixel(tiempo), y_pixel(valor))
        for tiempo, valor in zip(tiempos, valores_graficados)
    ]
    polilinea = " ".join(f"{x:.2f},{y:.2f}" for x, y in puntos)
    relleno = " ".join(
        [f"{x0:.2f},{y1:.2f}"]
        + [f"{x:.2f},{y:.2f}" for x, y in puntos]
        + [f"{x1:.2f},{y1:.2f}"]
    )
    elementos = [
        f'<svg viewBox="0 0 {ancho} {alto}" role="img" aria-label="{escape(titulo)}">',
        f'<rect x="0" y="0" width="{ancho}" height="{alto}" fill="#ffffff"/>',
        f'<text class="plot-title" x="{(x0+x1)/2}" y="25" text-anchor="middle">{escape(titulo)}</text>',
    ]

    for indice in range(5):
        proporcion = indice / 4
        x = x0 + proporcion * (x1 - x0)
        t = t_min + proporcion * (t_max - t_min)
        y = y1 - proporcion * (y1 - y0)
        valor = y_min + proporcion * (y_max - y_min)
        elementos.extend((
            f'<line class="grid" x1="{x:.2f}" y1="{y0}" x2="{x:.2f}" y2="{y1}"/>',
            f'<line class="grid" x1="{x0}" y1="{y:.2f}" x2="{x1}" y2="{y:.2f}"/>',
            f'<text class="tick" x="{x:.2f}" y="{y1+20}" text-anchor="middle">{t:.3g}</text>',
            f'<text class="tick" x="{x0-10}" y="{y+4:.2f}" text-anchor="end">{valor:.3g}</text>',
        ))

    elementos.extend((
        f'<line class="axis" x1="{x0}" y1="{y1}" x2="{x1}" y2="{y1}"/>',
        f'<line class="axis" x1="{x0}" y1="{y0}" x2="{x0}" y2="{y1}"/>',
        f'<polygon points="{relleno}" fill="{color}" fill-opacity="0.10"/>',
        f'<polyline points="{polilinea}" fill="none" stroke="{color}" stroke-width="2.6" stroke-linejoin="round" stroke-linecap="round"/>',
    ))
    for x, y in puntos:
        elementos.append(
            f'<circle cx="{x:.2f}" cy="{y:.2f}" r="2.2" fill="#ffffff" stroke="{color}" stroke-width="1.5"/>'
        )
    elementos.extend((
        f'<text class="axis-label" x="{(x0+x1)/2}" y="{alto-12}" text-anchor="middle">Tiempo t</text>',
        f'<text class="axis-label" transform="translate(18 {(y0+y1)/2}) rotate(-90)" text-anchor="middle">{escape(eje_y)}</text>',
    ))
    if estado == "divergente" and tiempo_fallo is not None:
        prefijo = "Escala recortada · " if recorte_escala else "Integración detenida · "
        etiqueta = f"{prefijo}t ≈ {tiempo_fallo:.4g}"
        elementos.extend((
            f'<rect x="{x1-236}" y="{y0+8}" width="226" height="28" rx="5" fill="#fff1f2" stroke="#fecdd3"/>',
            f'<text class="failure-label" x="{x1-123}" y="{y0+27}" text-anchor="middle">{escape(etiqueta)}</text>',
        ))
    elementos.append("</svg>")
    return "".join(elementos)

def crear_reporte_html(modelo, matriculas):
    _, N, P, S_base = modelo["resultados"]
    filas_base = "".join(
        "<tr>" + "".join(
            f"<td>{escape(str(valor))}</td>" for valor in (
                f"Integrante {clave[1:]}", "".join(map(str, matriculas[clave])),
                modelo["resultados"][0][clave], f"{N[clave]:.7f}",
                f"{P[clave]:.7f}", S_base[clave],
            )
        ) + "</tr>"
        for clave in ("M1", "M2", "M3", "M4")
    )
    variables = modelo["variables_a"]
    filas_a = "".join(
        f"<tr><td>{escape(clave)}</td><td>{escape(SIGNIFICADOS_A[clave])}</td>"
        f"<td><code>{escape(datos['formula'])}</code></td><td>{datos['valor']:.10g}</td></tr>"
        for clave, datos in variables.items()
    )
    b, c, d, e, f, g = (
        modelo["documento_b"], modelo["documento_c"], modelo["documento_d"],
        modelo["documento_e"], modelo["documento_f"], modelo["documento_g"],
    )
    resultados_documentos = (
        ("B", "Ar", SIGNIFICADOS_B["Ar"], b["Ar"]),
        ("B", "Br", SIGNIFICADOS_B["Br"], b["Br"]),
        ("B", "Cr", SIGNIFICADOS_B["Cr"], b["Cr"]),
        ("B", "Dr", SIGNIFICADOS_B["Dr"], b["Dr"]),
        ("C", "a_s", SIGNIFICADOS_C["a_s"], c["a_s"]),
        ("C", "b_s", SIGNIFICADOS_C["b_s"], c["b_s"]),
        ("C", "c_s", SIGNIFICADOS_C["c_s"], c["c_s"]),
        ("C", "d_s", SIGNIFICADOS_C["d_s"], c["d_s"]),
        ("D", "k", SIGNIFICADOS_D["k"], d["k"]),
        ("D", "a0", SIGNIFICADOS_D["a0"], d["a0"]),
        ("D", "m", SIGNIFICADOS_D["m"], d["m"]),
        ("D", "n", SIGNIFICADOS_D["n"], d["n"]),
        ("D", "p", SIGNIFICADOS_D["p"], d["p"]),
        ("D", "q", SIGNIFICADOS_D["q"], d["q"]),
        ("D", "L(3)", SIGNIFICADOS_D["L(3)"], d["L3"]),
        ("E", "r_sigma", SIGNIFICADOS_E["r_sigma"], e["r_sigma"]),
        ("E", "sigma", SIGNIFICADOS_E["sigma"], e["sigma"]),
        ("E", "tau", SIGNIFICADOS_E["tau"], e["tau"]),
        ("E", "kappa", SIGNIFICADOS_E["kappa"], e["kappa"]),
        ("E", "omega", SIGNIFICADOS_E["omega"], e["omega"]),
        ("E", "Y(0)", SIGNIFICADOS_E["Y(0)"], e["Y0"]),
        ("E", "Y(3)", SIGNIFICADOS_E["Y(3)"], e["Y_final"]),
        ("F", "chi", SIGNIFICADOS_F["chi"], f["chi"]),
        ("F", "psi", SIGNIFICADOS_F["psi"], f["psi"]),
        ("F", "xi", SIGNIFICADOS_F["xi"], f["xi"]),
        ("F", "Estado de integración", SIGNIFICADOS_F["estado"], f["estado"]),
        ("F", "Menor aporte Y/xi", SIGNIFICADOS_F["fuente_minima"], f["fuente_minima"]),
        ("F", "Mayor aporte Y/xi", SIGNIFICADOS_F["Y/xi"], f["fuente_maxima"]),
        ("F", "Capacidad máxima", SIGNIFICADOS_F["capacidad"], f["capacidad_maxima"]),
        ("F", "Temperatura inicial θ(0)", SIGNIFICADOS_F["theta(0)"], f["theta0"]),
        ("F", "Temperatura final θ(3)", SIGNIFICADOS_F["theta(3)"], f["theta_final"]),
        ("F", "Tiempo detenido", "Tiempo aproximado donde RK4 dejó de avanzar.", f["tiempo_fallo"]),
        ("F", "Diagnóstico", "Motivo por el que se detuvo la integración.", f["causa"] or "Ninguno."),
    )
    if g["estado"] == "completo":
        resultados_documentos += (
            ("G", "Integral de Y", SIGNIFICADOS_G["Integral de Y"], g["integral_Y"]),
            ("G", "Temperatura θ(1)", SIGNIFICADOS_G["theta(1)"], g["theta_1"]),
            ("G", "Nivel operativo U*", SIGNIFICADOS_G["U*"], g["U_star"]),
            ("G", "Índice general S", SIGNIFICADOS_G["S"], g["S"]),
        )
    else:
        resultados_documentos += (
            ("G", "Estado", g["causa"], "No disponible"),
        )
    filas_modelos_lista = []
    for documento, nombre, descripcion, valor in resultados_documentos:
        if documento == "F" and nombre == "Estado de integración":
            estado_html = (
                '<span class="badge badge-danger">Divergente</span>'
                if valor == "divergente"
                else '<span class="badge badge-success">Completo</span>'
            )
        elif documento == "G" and nombre == "Estado":
            estado_html = '<span class="badge badge-warning">No disponible</span>'
            valor = estado_html
        else:
            estado_html = ""
        if estado_html and documento == "F":
            valor_html = estado_html
        elif estado_html and documento == "G":
            valor_html = valor
        elif isinstance(valor, str):
            valor_html = escape(valor)
        elif valor is None:
            valor_html = '<span class="muted">No disponible</span>'
        else:
            valor_html = f"{valor:.10g}"
        filas_modelos_lista.append(
            f"<tr><td>Documento {escape(documento)}</td><td>{escape(nombre)}</td>"
            f"<td>{escape(str(descripcion))}</td><td>{valor_html}</td></tr>"
        )
    filas_modelos = "".join(filas_modelos_lista)
    desglose_cs = c["desglose_c_s"]
    filas_cs_tabla = "".join(
        f"<tr><td>{tiempo:g}</td><td>{valor:.4g}</td></tr>"
        for tiempo, valor in zip(desglose_cs["nodos"], desglose_cs["valores"])
    )
    filas_cs_desglose = "".join(
        f"<tr><td>{escape(etiqueta)}</td><td>{valor}</td></tr>"
        for etiqueta, valor in (
            ("Número de intervalos n", str(desglose_cs["num_intervalos"])),
            ("Paso h", f"{desglose_cs['h']:g}"),
            ("Σ extremos (l0 + l4)", f"{desglose_cs['suma_extremos']:.10g}"),
            ("Σ impares (l1 + l3)", f"{desglose_cs['suma_impares']:.10g}"),
            ("Σ pares internos (l2)", f"{desglose_cs['suma_pares']:.10g}"),
            (
                "Σ ponderada = extremos + 4·impares + 2·pares",
                f"{desglose_cs['suma_ponderada']:.10g}",
            ),
            ("c_s = (h/3)·Σ ponderada", f"{desglose_cs['resultado']:.10g}"),
        )
    )
    estado_f_texto = "Divergente" if f["estado"] == "divergente" else "Completo"
    estado_f_clase = "badge-danger" if f["estado"] == "divergente" else "badge-success"
    estado_g_texto = "No disponible" if g["estado"] != "completo" else "Calculado"
    estado_g_clase = "badge-warning" if g["estado"] != "completo" else "badge-success"
    equipo_texto = ", ".join(
        "".join(map(str, matriculas[clave])) for clave in ("M1", "M2", "M3", "M4")
    )
    filas_sensibilidad = []
    for dato in modelo["sensibilidad"]["pasos"]:
        diferencia = "—" if dato["diferencia_S"] is None else f"{dato['diferencia_S']:.3g}"
        theta_final = "No disponible" if dato["theta_final"] is None else f"{dato['theta_final']:.8g}"
        indice = "No disponible" if dato["S"] is None else f"{dato['S']:.8g}"
        fallo = "—" if dato["tiempo_fallo_F"] is None else f"{dato['tiempo_fallo_F']:.5g}"
        filas_sensibilidad.append(
            f"<tr><td>{dato['h']:.3g}</td><td>{dato['Y_final']:.8g}</td>"
            f"<td>{theta_final}</td><td>{indice}</td><td>{diferencia}</td>"
            f"<td>{escape(dato['estado_F'])}</td><td>{fallo}</td></tr>"
        )
    filas_sensibilidad = "".join(filas_sensibilidad)
    filas_tolerancias = "".join(
        f"<tr><td>{dato['tolerancia']:.0e}</td><td>{dato['C']:.10g}</td>"
        f"<td>{dato['iteraciones']}</td><td>{abs(dato['residuo']):.3g}</td>"
        f"<td>{'No disponible' if dato['S'] is None else format(dato['S'], '.10g')}</td>"
        f"<td>{escape(dato['estado_F'])}</td></tr>"
        for dato in modelo["sensibilidad"]["tolerancias"]
    )
    residuo_n = d["n"]**3 - 7*d["n"] + 6
    residuo_u = None
    if g["estado"] == "completo":
        U_star = g["U_star"]
        residuo_u = c["a_s"]*U_star**3 - c["b_s"]*U_star**2 + c["c_s"]*U_star - c["d_s"]
    C = variables["C"]["valor"]
    residuo_c = C**3 - (P["M1"] + P["M2"])*C**2 + S_base["M3"]*C - 5
    codigo_fuente = escape(Path(__file__).read_text(encoding="utf-8"))
    return f"""<!doctype html>
<html lang="es"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>Informe técnico PIA</title>
<style>
@page{{size:A4;margin:16mm 15mm}}
:root{{color-scheme:light;--ink:#172b3b;--muted:#607487;--navy:#102a43;--cyan:#087e8b;--line:#d5e0e8;--paper:#f4f7fa;--white:#fff}}
*{{box-sizing:border-box}}
body{{margin:0;background:#e8eef3;color:var(--ink);font:14px/1.58 "Segoe UI",Arial,sans-serif}}
.report{{max-width:1040px;margin:32px auto;padding:0 28px 44px}}
.report-head{{position:relative;overflow:hidden;padding:38px 42px 32px;background:var(--navy);color:#fff;border-radius:12px 12px 4px 4px}}
.report-head:after{{content:"";position:absolute;right:-68px;top:-112px;width:300px;height:300px;border:1px solid #ffffff25;border-radius:50%;box-shadow:0 0 0 34px #ffffff0d,0 0 0 68px #ffffff08}}
.brand{{display:flex;align-items:center;gap:12px;color:#a8e8ec;font-size:12px;font-weight:700;letter-spacing:.04em}}
.brand-mark{{width:28px;height:28px;display:grid;place-items:center;border:1px solid #53c5ce;border-radius:7px;color:#78e1e6;font-size:11px}}
.report-head h1{{position:relative;z-index:1;margin:26px 0 8px;max-width:700px;color:#fff;font-size:30px;line-height:1.17;letter-spacing:-.025em}}
.report-head p{{position:relative;z-index:1;max-width:700px;margin:0;color:#c5d3df;font-size:14px}}
.team-line{{position:relative;z-index:1;margin-top:24px;padding-top:16px;border-top:1px solid #ffffff2a;color:#c5d3df;font-size:12px}}
.snapshot{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:12px;margin:18px 0 26px}}
.snapshot-card{{min-height:92px;padding:16px 18px;background:var(--white);border:1px solid var(--line);border-radius:8px}}
.snapshot-label{{display:block;margin-bottom:8px;color:var(--muted);font-size:11px;font-weight:700}}
.snapshot-value{{font-size:19px;font-weight:700;letter-spacing:-.02em;color:var(--navy)}}
.report-section{{margin:16px 0;padding:24px 26px;background:var(--white);border:1px solid var(--line);border-radius:8px}}
h2{{margin:0 0 8px;color:var(--navy);font-size:19px;line-height:1.3}}
h3{{margin:22px 0 8px;color:var(--navy);font-size:15px}}
h2,h3{{break-after:avoid;page-break-after:avoid}}
p{{margin:8px 0 14px;color:#40566a}}
.section-rule{{width:42px;height:3px;margin:0 0 16px;background:var(--cyan);border-radius:2px}}
table{{width:100%;margin:14px 0 4px;border-collapse:separate;border-spacing:0;border:1px solid var(--line);border-radius:6px;overflow:hidden;font-size:12px}}
thead{{display:table-header-group}}
th,td{{padding:9px 10px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top}}
th{{background:#eaf1f5;color:#264257;font-size:11px;font-weight:700}}
tbody tr:nth-child(even){{background:#f7fafc}}
tbody tr:last-child td{{border-bottom:0}}
tr{{break-inside:avoid;page-break-inside:avoid}}
code{{color:#075b67;font-family:Consolas,"Courier New",monospace;font-size:.92em;overflow-wrap:anywhere}}
.badge{{display:inline-flex;align-items:center;padding:3px 9px;border-radius:20px;font-size:11px;font-weight:700;white-space:nowrap}}
.badge-danger{{background:#fff0ee;color:#aa3024;border:1px solid #ffd1ca}}
.badge-warning{{background:#fff6e3;color:#8a5a08;border:1px solid #f4dda9}}
.badge-success{{background:#eaf8f1;color:#176b46;border:1px solid #bfe7d2}}
.muted{{color:#75879a;font-style:italic}}
.chart-card{{margin:14px 0 18px;padding:10px 12px 4px;border:1px solid var(--line);border-radius:7px;background:#fff;break-inside:avoid;page-break-inside:avoid}}
svg{{display:block;width:100%;height:auto}}
svg .plot-title{{fill:#203b50;font:600 14px "Segoe UI",Arial,sans-serif}}
svg .grid{{stroke:#dbe5ed;stroke-width:1;stroke-dasharray:3 5}}
svg .axis{{stroke:#71869a;stroke-width:1.2}}
svg .tick{{fill:#63798d;font:11px "Segoe UI",Arial,sans-serif}}
svg .axis-label{{fill:#344e62;font:12px "Segoe UI",Arial,sans-serif}}
svg .failure-label{{fill:#a52621;font:700 11px "Segoe UI",Arial,sans-serif}}
.method-list{{margin:10px 0 0;padding-left:20px;color:#40566a}}
.method-list li{{margin:7px 0}}
.source-section{{break-before:page;page-break-before:always}}
pre{{margin:12px 0 0;padding:16px;background:#f1f5f8;border:1px solid var(--line);border-radius:6px;white-space:pre-wrap;overflow-wrap:anywhere;font-size:7.5pt;line-height:1.35}}
.report-foot{{padding:14px 2px;color:var(--muted);font-size:11px}}
@media(max-width:720px){{.report{{padding:0 12px 28px;margin:16px auto}}.report-head{{padding:28px 24px}}.report-head h1{{font-size:25px}}.report-section{{padding:18px 16px;overflow-x:auto}}.snapshot{{grid-template-columns:1fr}}table{{min-width:620px}}}}
@media print{{body{{background:#fff;font-size:10pt}}.report{{max-width:none;margin:0;padding:0}}.report-head,.badge,.snapshot-card,th,.chart-card{{-webkit-print-color-adjust:exact;print-color-adjust:exact}}.report-head{{border-radius:0}}.snapshot-card,.report-section,.chart-card{{box-shadow:none}}.report-section{{margin:10px 0;padding:18px}}.source-section{{break-before:page}}a{{color:inherit;text-decoration:none}}}}
</style></head>
<body><main class="report">
<header class="report-head">
<div class="brand"><span class="brand-mark">MN</span><span>PIA · Métodos Numéricos</span></div>
<h1>Análisis de respuesta y estabilidad de una plataforma distribuida</h1>
<p>Informe técnico del Proyecto Integrador de Aprendizaje. Parámetros reconstruidos, métodos aplicados y límites observados en la simulación.</p>
<div class="team-line">Equipo: {equipo_texto}</div>
</header>
<div class="snapshot">
<div class="snapshot-card"><span class="snapshot-label">Equipo analizado</span><span class="snapshot-value">4 matrículas</span></div>
<div class="snapshot-card"><span class="snapshot-label">Documento F · integración térmica</span><span class="badge {estado_f_clase}">{estado_f_texto}</span></div>
<div class="snapshot-card"><span class="snapshot-label">Documento G · índice global</span><span class="badge {estado_g_clase}">{estado_g_texto}</span></div>
</div>
<section class="report-section">
<h2>Datos del equipo</h2><div class="section-rule"></div>
<p><b>M</b> es la matrícula escrita como número; <b>N</b> es su escala normalizada; <b>P</b> es el promedio de dígitos y <b>S base</b> es su suma ponderada. S base no es el índice del Documento G.</p>
<table><thead><tr><th>Integrante</th><th>Matrícula</th><th>M</th><th>N</th><th>P</th><th>S base</th></tr></thead><tbody>{filas_base}</tbody></table>
</section>
<section class="report-section">
<h2>Documento A · variables operativas</h2><div class="section-rule"></div>
<p>Las descripciones explican la función de cada resultado; la operación técnica conserva la relación usada por el modelo.</p>
<table><thead><tr><th>Clave</th><th>Interpretación</th><th>Operación técnica</th><th>Resultado</th></tr></thead><tbody>{filas_a}</tbody></table>
</section>
<section class="report-section">
<h2>Documentos B–G · parámetros reconstruidos</h2><div class="section-rule"></div>
<table><thead><tr><th>Documento</th><th>Parámetro</th><th>Interpretación</th><th>Resultado</th></tr></thead><tbody>{filas_modelos}</tbody></table>
</section>
<section class="report-section">
<h2>Documento C · desglose de c_s con Simpson 1/3 compuesto</h2><div class="section-rule"></div>
<p>La tabla l(t) tiene 5 puntos y 4 intervalos, por lo que Simpson 1/3 compuesto es aplicable (n par). Pesos: 1, 4, 2, 4, 1 con paso h = {desglose_cs['h']:g}.</p>
<table><thead><tr><th>t</th><th>l(t)</th></tr></thead><tbody>{filas_cs_tabla}</tbody></table>
<table><thead><tr><th>Concepto</th><th>Valor</th></tr></thead><tbody>{filas_cs_desglose}</tbody></table>
</section>
<section class="report-section">
<h2>Trayectorias de respuesta</h2><div class="section-rule"></div>
<p>Los ejes identifican el tiempo y la magnitud representada. Los puntos indican los pasos de integración. En F, el tramo extremo se recorta solo en la gráfica para conservar visibles los cambios previos; los valores del modelo no se modifican.</p>
<div class="chart-card">{crear_svg_trayectoria(e["historial"], "Documento E: respuesta del sistema", "#087e8b", "Respuesta Y(t)")}</div>
<div class="chart-card">{crear_svg_trayectoria(f["historial"], "Documento F: temperatura del sistema", "#bd453a", "Temperatura θ(t)", f["estado"], f["tiempo_fallo"])}</div>
</section>
<section class="report-section">
<h2>Sensibilidad numérica</h2><div class="section-rule"></div>
<p>La tabla compara tres tamaños de paso RK4 y tres tolerancias para Newton-Raphson en C. La temperatura final y G se marcan no disponibles cuando F no completa una trayectoria finita.</p>
<table><thead><tr><th>Paso h</th><th>Y(3)</th><th>θ(3)</th><th>Índice G</th><th>Diferencia</th><th>Estado F</th><th>Fallo t</th></tr></thead><tbody>{filas_sensibilidad}</tbody></table>
<h3>Tolerancia de Newton-Raphson para C</h3>
<table><thead><tr><th>Tolerancia</th><th>C</th><th>Iteraciones</th><th>Residuo</th><th>Índice G</th><th>Estado F</th></tr></thead><tbody>{filas_tolerancias}</tbody></table>
</section>
<section class="report-section">
<h2>Métodos, criterios y supuestos</h2><div class="section-rule"></div>
<ul class="method-list">
<li>Los métodos numéricos se implementan en el código: bisección y aislamiento de raíces polinómicas, Newton-Raphson para C, interpolación de Lagrange y Newton hacia adelante, mínimos cuadrados, Gauss-Jordan con pivoteo, Simpson 1/3 compuesto (sobre funciones y sobre tablas, con desglose de sumatorias), Simpson 3/8 y Runge-Kutta de cuarto orden.</li>
<li>Para el modelo E se selecciona r_sigma = 1.2 entre las raíces admisibles 1.2 y 4.0; se usa la raíz menor para adoptar el decaimiento más conservador.</li>
<li>Se conserva la asignación del expediente: kappa = L(3), donde L se forma con los coeficientes del Documento D; omega es la integral de la tabla v(t) en [0,3], calculada con Simpson 3/8.</li>
<li>En F se propone theta(0) = 25 °C como temperatura ambiente/de reposo previa a la carga de trabajo, conforme a la justificación del equipo.</li>
<li>Se conserva sin normalización la ecuación original theta' = chi*theta*(psi-theta) - Y/xi. Si Y/xi supera chi*psi²/4 durante todo el intervalo, el aporte con signo negativo domina al término logístico y theta' es negativa; cuando theta cae por debajo de cero, el término logístico también es negativo y la magnitud crece hasta que RK4 deja de obtener valores finitos. Se reporta el diagnóstico, sin inventar theta(3), theta(1) ni G.</li>
<li>Para U* se resuelve la ecuación cúbica de saturación y se elige la menor raíz real positiva. El criterio prioriza el nivel operativo positivo más bajo.</li>
<li>El análisis de sensibilidad compara pasos de RK4 y tolerancias de Newton. Si el modelo no produce una trayectoria finita, el índice global no se define para esa instancia y se marca como no disponible.</li>
</ul></section>
<section class="report-section">
<h2>Interpretación y validación</h2><div class="section-rule"></div>
<p>El índice G solo se calcula cuando θ(1) está disponible. Si F diverge antes de t=1, el resultado correcto es informar que el índice no puede calcularse con el modelo y los datos originales; no se sustituye por un resultado basado en una ecuación modificada.</p>
<table><thead><tr><th>Comprobación</th><th>Resultado</th></tr></thead><tbody>
<tr><td>Pendiente k con la tabla del documento D</td><td>{d["k"]:.10g} (referencia manual: 1.71)</td></tr>
<tr><td>Integral m, polinomio integrado exactamente por Simpson</td><td>{d["m"]:.10g} (valor analítico: 15)</td></tr>
<tr><td>Residuo de la raíz n de D</td><td>{residuo_n:.3g}</td></tr>
<tr><td>Acumulado c_s con Simpson 1/3 compuesto</td><td>{c["c_s"]:.10g} (n={desglose_cs["num_intervalos"]}, h={desglose_cs["h"]:g}, Σ ponderada={desglose_cs["suma_ponderada"]:.10g})</td></tr>
<tr><td>Residuo de Newton-Raphson para C</td><td>{residuo_c:.3g}</td></tr>
<tr><td>Residuo de la raíz U* de saturación</td><td>{"No disponible: Documento F/G no computable" if residuo_u is None else f"{residuo_u:.3g}"}</td></tr>
</tbody></table>
<p>Los documentos A–F y los cálculos intermedios se reconstruyen desde las matrículas. La evaluación de G queda condicionada a que la ecuación F produzca θ(1) finita; cualquier divergencia se informa sin modificar el modelo original.</p>
<p><b>Reproducibilidad.</b> Programa PIA.py; requiere Python 3.10 o posterior y Tkinter, sin dependencias externas. Ejecutar con <code>python PIA.py</code>, capturar cuatro matrículas y usar Exportar informe. El HTML se puede imprimir a PDF.</p>
</section>
<section class="report-section source-section">
<h2>Código fuente completo</h2><div class="section-rule"></div><p>Versión ejecutable usada para generar este informe.</p>
<pre><code>{codigo_fuente}</code></pre>
</section>
<footer class="report-foot">Proyecto Integrador de Aprendizaje · Métodos Numéricos · Informe técnico</footer>
</main></body></html>"""

class Aplicacion(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("PIA · Métodos Numéricos")
        self.geometry("1240x800")
        self.minsize(1000, 720)
        self.colores = {
            "bg_principal": "#0b0f19",
            "bg_tarjeta": "#111827",
            "sidebar": "#080d17",
            "texto_oscuro": "#eef3f8",
            "texto_secundario": "#9aa9b9",
            "primario": "#38bdf8",
            "primario_hover": "#7dd3fc",
            "acento": "#38bdf8",
            "barra_superior": "#080d17",
            "borde": "#374151",
            "linea_tabla": "#1f2937",
            "grid": "#263244",
        }
        self.configure(bg=self.colores["bg_principal"])
        self.matriculas = {}
        self.resultados = None
        self.modelo = None
        self.entradas = []
        self.active_view = "capture"
        self._nav_buttons = {}
        self._scroll_canvas = None
        self._scroll_window = None
        self._build_layout()
        self._build_sidebar()
        self.bind_all("<MouseWheel>", self._on_mousewheel, add="+")
        self.bind_all("<Button-4>", lambda event: self._on_mousewheel(event, -1), add="+")
        self.bind_all("<Button-5>", lambda event: self._on_mousewheel(event, 1), add="+")
        self.mostrar_captura()

    def _build_layout(self):
        estilo = ttk.Style(self)
        if "clam" in estilo.theme_names():
            estilo.theme_use("clam")
        estilo.configure(
            "PIA.Vertical.TScrollbar",
            background="#263244", troughcolor="#0b0f19",
            bordercolor="#111827", arrowcolor="#9aa9b9",
            lightcolor="#263244", darkcolor="#263244", width=10,
        )
        estilo.map(
            "PIA.Vertical.TScrollbar",
            background=[("active", "#38bdf8")],
        )
        estilo.configure(
            "PIA.Horizontal.TScrollbar",
            background="#263244", troughcolor="#0b0f19",
            bordercolor="#111827", arrowcolor="#9aa9b9",
            lightcolor="#263244", darkcolor="#263244", width=10,
        )
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)
        shell = tk.Frame(self, bg=self.colores["bg_principal"])
        shell.grid(row=0, column=0, sticky="nsew")
        shell.grid_rowconfigure(0, weight=1)
        shell.grid_columnconfigure(1, weight=1)

        self.sidebar = tk.Frame(
            shell, bg=self.colores["sidebar"], width=252,
            highlightthickness=1, highlightbackground=self.colores["borde"],
        )
        self.sidebar.grid(row=0, column=0, sticky="ns")
        self.sidebar.grid_propagate(False)

        workspace = tk.Frame(shell, bg=self.colores["bg_principal"])
        workspace.grid(row=0, column=1, sticky="nsew")
        workspace.grid_rowconfigure(0, weight=1)
        workspace.grid_columnconfigure(0, weight=1)
        self.content_canvas = tk.Canvas(
            workspace, bg=self.colores["bg_principal"], highlightthickness=0,
        )
        self.content_scrollbar = ttk.Scrollbar(
            workspace, orient="vertical", command=self.content_canvas.yview,
            style="PIA.Vertical.TScrollbar",
        )
        self.content_canvas.grid(row=0, column=0, sticky="nsew")
        self.content_scrollbar.grid(row=0, column=1, sticky="ns")
        self.content_canvas.configure(yscrollcommand=self.content_scrollbar.set)
        self.content = tk.Frame(self.content_canvas, bg=self.colores["bg_principal"])
        self._scroll_window = self.content_canvas.create_window(
            (0, 0), window=self.content, anchor="nw",
        )
        self.content_canvas._pia_scroll_canvas = self.content_canvas
        self.content.bind(
            "<Configure>",
            lambda event: self.content_canvas.configure(
                scrollregion=self.content_canvas.bbox("all")
            ),
        )
        self.content_canvas.bind(
            "<Configure>",
            lambda event: self.content_canvas.itemconfigure(
                self._scroll_window, width=event.width
            ),
        )

    def _build_sidebar(self):
        brand = tk.Frame(self.sidebar, bg=self.colores["sidebar"], padx=20, pady=22)
        brand.pack(fill="x")
        tk.Label(
            brand, text="PIA · MÉTODOS NUMÉRICOS",
            bg=self.colores["sidebar"], fg=self.colores["acento"],
            font=("Segoe UI", 10, "bold"), anchor="w",
        ).pack(fill="x")
        tk.Label(
            brand, text="Análisis de sistemas\ndistribuidos",
            bg=self.colores["sidebar"], fg=self.colores["texto_oscuro"],
            font=("Segoe UI", 15, "bold"), justify="left", anchor="w",
        ).pack(fill="x", pady=(9, 0))
        tk.Frame(brand, bg=self.colores["borde"], height=1).pack(fill="x", pady=(20, 0))

        navigation = tk.Frame(self.sidebar, bg=self.colores["sidebar"], padx=11)
        navigation.pack(fill="both", expand=True, pady=(4, 8))
        self._sidebar_section(navigation, "VISTA GENERAL")
        rutas = (
            ("results", "Resultados base", self.mostrar_resultados),
            ("analysis", "Análisis y gráficas", self.mostrar_analisis),
            ("operations", "Operaciones detalladas", self.mostrar_operaciones),
        )
        for key, label, command in rutas:
            self._add_nav_button(navigation, key, label, command)

        self._sidebar_section(navigation, "DOCUMENTOS")
        for key, label, command in (
            ("doc_a", "A · Estado operativo", self.mostrar_documento_a),
            ("doc_b", "B · Respuesta bajo carga", self.mostrar_documento_b),
            ("doc_c", "C · Saturación", self.mostrar_documento_c),
            ("doc_d", "D · Balance acumulado", self.mostrar_documento_d),
            ("doc_e", "E · Rendimiento", self.mostrar_documento_e),
            ("doc_f", "F · Temperatura", self.mostrar_documento_f),
            ("doc_g", "G · Índice global", self.mostrar_documento_g),
        ):
            self._add_nav_button(navigation, key, label, command)

        actions = tk.Frame(
            self.sidebar, bg=self.colores["sidebar"], padx=16, pady=16,
        )
        actions.pack(side="bottom", fill="x")
        tk.Frame(actions, bg=self.colores["borde"], height=1).pack(fill="x", pady=(0, 14))
        tk.Button(
            actions, text="Exportar informe HTML",
            command=self.exportar_reporte, anchor="w",
            bg=self.colores["acento"], fg="#06202c",
            activebackground=self.colores["primario_hover"], activeforeground="#06202c",
            relief="flat", padx=13, pady=11, font=("Segoe UI", 10, "bold"),
            cursor="hand2",
        ).pack(fill="x")
        tk.Button(
            actions, text="Nueva captura", command=self.mostrar_captura,
            anchor="w", bg=self.colores["sidebar"], fg=self.colores["texto_oscuro"],
            activebackground="#172235", activeforeground="white",
            relief="flat", padx=13, pady=10, font=("Segoe UI", 10),
            cursor="hand2",
        ).pack(fill="x", pady=(7, 0))

    def _sidebar_section(self, parent, title):
        tk.Label(
            parent, text=title, bg=self.colores["sidebar"],
            fg=self.colores["texto_secundario"], anchor="w",
            font=("Segoe UI", 8, "bold"),
        ).pack(fill="x", padx=11, pady=(11, 5))

    def _add_nav_button(self, parent, key, label, command):
        row = tk.Frame(parent, bg=self.colores["sidebar"])
        row.pack(fill="x", pady=1)
        indicator = tk.Frame(row, bg=self.colores["sidebar"], width=3)
        indicator.pack(side="left", fill="y")
        button = tk.Button(
            row, text=label, command=command, anchor="w",
            bg=self.colores["sidebar"], fg=self.colores["texto_secundario"],
            activebackground="#172235", activeforeground=self.colores["texto_oscuro"],
            relief="flat", padx=11, pady=8,
            font=("Segoe UI", 9), cursor="hand2",
        )
        button.pack(side="left", fill="x", expand=True)
        self._nav_buttons[key] = (indicator, button)

    def _set_active(self, key):
        self.active_view = key
        for view_key, (indicator, button) in self._nav_buttons.items():
            active = view_key == key
            indicator.configure(bg=self.colores["acento"] if active else self.colores["sidebar"])
            button.configure(
                bg="#132235" if active else self.colores["sidebar"],
                fg=self.colores["texto_oscuro"] if active else self.colores["texto_secundario"],
                font=("Segoe UI", 9, "bold" if active else "normal"),
            )

    def _on_mousewheel(self, event, direction=None):
        widget = self.winfo_containing(event.x_root, event.y_root)
        while widget is not None:
            canvas = getattr(widget, "_pia_scroll_canvas", None)
            if canvas is not None and canvas.winfo_exists():
                delta = direction if direction is not None else (-1 if event.delta > 0 else 1)
                canvas.yview_scroll(delta, "units")
                return "break"
            widget = getattr(widget, "master", None)
        return None

    def limpiar_pantalla(self):
        for widget in self.content.winfo_children():
            widget.destroy()
        self.content_canvas.yview_moveto(0)

    def crear_encabezado(self, titulo, subtitulo=""):
        header_frame = tk.Frame(self.content, bg=self.colores["bg_principal"])
        header_frame.pack(fill="x", padx=34, pady=(28, 16))
        tk.Frame(header_frame, bg=self.colores["acento"], width=3).pack(
            side="left", fill="y", padx=(0, 15),
        )
        title_block = tk.Frame(header_frame, bg=self.colores["bg_principal"])
        title_block.pack(side="left", fill="x", expand=True)
        tk.Label(
            title_block, text=titulo,
            font=("Segoe UI", 23, "bold"), fg=self.colores["texto_oscuro"],
            bg=self.colores["bg_principal"], anchor="w",
        ).pack(anchor="w")
        if subtitulo:
            tk.Label(
                title_block, text=subtitulo, font=("Segoe UI", 10),
                fg=self.colores["texto_secundario"], bg=self.colores["bg_principal"],
                anchor="w", justify="left", wraplength=900,
            ).pack(anchor="w", pady=(5, 0))

    def agregar_navegacion(self):
        self._set_active(self.active_view)

    def crear_tarjetas_kpi(self, elementos):
        contenedor = tk.Frame(self.content, bg=self.colores["bg_principal"])
        contenedor.pack(fill="x", padx=34, pady=(0, 18))
        for columna in range(len(elementos)):
            contenedor.grid_columnconfigure(columna, weight=1, uniform="kpi")
        for columna, (titulo, valor, detalle, tono) in enumerate(elementos):
            tarjeta = tk.Frame(
                contenedor, bg=self.colores["bg_tarjeta"],
                highlightthickness=1, highlightbackground=self.colores["borde"],
                padx=16, pady=13,
            )
            tarjeta.grid(row=0, column=columna, sticky="nsew", padx=(0 if columna == 0 else 8, 0))
            tk.Label(
                tarjeta, text=titulo, bg=self.colores["bg_tarjeta"],
                fg=self.colores["texto_secundario"], font=("Segoe UI", 9),
                anchor="w",
            ).pack(fill="x")
            tk.Label(
                tarjeta, text=valor, bg=self.colores["bg_tarjeta"],
                fg=tono, font=("Segoe UI", 18, "bold"), anchor="w",
            ).pack(fill="x", pady=(5, 2))
            tk.Label(
                tarjeta, text=detalle, bg=self.colores["bg_tarjeta"],
                fg=self.colores["texto_secundario"], font=("Segoe UI", 8),
                anchor="w", justify="left", wraplength=260,
            ).pack(fill="x")

    def mostrar_captura(self):
        self.limpiar_pantalla()
        self._set_active("capture")
        self.crear_encabezado(
            "Iniciar análisis",
            "Capture las matrículas del equipo para reconstruir los parámetros del expediente.",
        )
        marco_central = tk.Frame(
            self.content, bg=self.colores["bg_tarjeta"],
            highlightthickness=1, highlightbackground=self.colores["borde"],
        )
        marco_central.pack(fill="x", padx=34, pady=(15, 30))
        inner_frame = tk.Frame(marco_central, bg=self.colores["bg_tarjeta"], padx=32, pady=30)
        inner_frame.pack(fill="both", expand=True)

        tk.Label(
            inner_frame, text="PIA · MÉTODOS NUMÉRICOS",
            font=("Segoe UI", 9, "bold"), fg=self.colores["primario"],
            bg=self.colores["bg_tarjeta"]
        ).pack(anchor="w")

        tk.Label(
            inner_frame, text="Matrículas del equipo",
            font=("Segoe UI", 19, "bold"), fg=self.colores["texto_oscuro"],
            bg=self.colores["bg_tarjeta"]
        ).pack(anchor="w", pady=(6, 8))

        tk.Label(
            inner_frame, text="Use cuatro matrículas distintas de siete dígitos. El sistema las ordenará antes de asignar M1–M4.",
            font=("Segoe UI", 10), fg=self.colores["texto_secundario"],
            bg=self.colores["bg_tarjeta"], wraplength=700, justify="left"
        ).pack(anchor="w", pady=(0, 18))

        self.entradas = []
        self.entrada_labels = []
        for i in range(1, 5):
            fila_mat = tk.Frame(inner_frame, bg=self.colores["bg_tarjeta"])
            fila_mat.pack(fill="x", pady=5)
            lbl = tk.Label(
                fila_mat, text=f"M{i}", width=5, anchor="w",
                font=("Segoe UI", 10, "bold"), fg=self.colores["texto_secundario"],
                bg=self.colores["bg_tarjeta"]
            )
            lbl.pack(side="left")

            ent = tk.Entry(
                fila_mat, font=("Consolas", 13), width=24,
                justify="left", relief="flat", bd=0, borderwidth=0,
                bg="#0b1220", fg=self.colores["texto_oscuro"], insertbackground="white",
                highlightthickness=1, highlightbackground=self.colores["borde"],
                highlightcolor=self.colores["acento"],
            )
            ent.pack(side="left", fill="x", expand=True, ipady=8, padx=(8, 0))
            ent.bind("<FocusIn>", lambda event, entry=ent: entry.configure(
                highlightbackground=self.colores["acento"],
            ))
            ent.bind("<FocusOut>", lambda event, entry=ent: self._validate_entry_visual(entry, True))
            ent.bind("<KeyRelease>", lambda event, entry=ent: self._validate_entry_visual(entry))
            ent.bind("<Return>", lambda event, index=i: self._advance_entry(index))
            self.entradas.append(ent)
            self.entrada_labels.append(lbl)

        btn_calcular = tk.Button(
            inner_frame, text="Calcular modelos", command=self.procesar,
            bg=self.colores["primario"], fg="white",
            activebackground=self.colores["primario_hover"], activeforeground="white",
            relief="flat", font=("Segoe UI", 10, "bold"), cursor="hand2",
            padx=18, pady=11
        )
        btn_calcular.pack(anchor="w", pady=(20, 0))
        tk.Label(
            inner_frame, text="Validación visual: siete dígitos numéricos por campo.",
            bg=self.colores["bg_tarjeta"], fg=self.colores["texto_secundario"],
            font=("Segoe UI", 8),
        ).pack(anchor="w", pady=(10, 0))

    def _advance_entry(self, index):
        if index < len(self.entradas):
            self.entradas[index].focus_set()
        else:
            self.procesar()
        return "break"

    def _validate_entry_visual(self, entry, final=False):
        valor = entry.get().strip()
        en_progreso = valor.isdigit() and len(valor) < 7
        valido = (
            not valor or (valor.isdigit() and len(valor) == 7)
            or (en_progreso and not final)
        )
        color = self.colores["borde"] if valido else "#f87171"
        entry.configure(highlightbackground=color)
        index = self.entradas.index(entry) if entry in self.entradas else None
        if index is not None:
            self.entrada_labels[index].configure(
                fg=self.colores["texto_secundario"] if valido else "#f87171"
            )

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

        valores_ordenados = sorted(valores)

        self.matriculas = {
            f"M{i}": [int(d) for d in val]
            for i, val in enumerate(valores_ordenados, start=1)
        }
        self.modelo = None
        self.resultados = None
        try:
            self.modelo = calcular_modelo(self.matriculas)
        except (ValueError, ZeroDivisionError, OverflowError) as error:
            self.modelo = None
            messagebox.showerror("No se pudo completar el modelo", str(error))
            return
        self.resultados = self.modelo["resultados"]
        self.mostrar_resultados()

    def calcular_documento_a(self):
        if not self.resultados:
            return {}
        return self.modelo["variables_a"]

    def crear_tabla_resultados(self, contenedor, filas_datos, encabezados):
        canvas_tabla = tk.Canvas(
            contenedor, bg=self.colores["bg_tarjeta"], highlightthickness=1,
            highlightbackground=self.colores["borde"],
        )
        canvas_tabla._pia_scroll_canvas = canvas_tabla
        scrollbar_vertical = ttk.Scrollbar(
            contenedor, orient="vertical", command=canvas_tabla.yview,
            style="PIA.Vertical.TScrollbar",
        )
        scrollbar_horizontal = ttk.Scrollbar(
            contenedor, orient="horizontal", command=canvas_tabla.xview,
            style="PIA.Horizontal.TScrollbar",
        )
        scrollable_frame = tk.Frame(canvas_tabla, bg=self.colores["bg_tarjeta"])
        fuente_celda = tkfont.Font(font=("Segoe UI", 10))
        fuente_numero = tkfont.Font(font=("Consolas", 10))
        fuente_encabezado = tkfont.Font(font=("Segoe UI", 9, "bold"))
        anchos = []
        for indice, encabezado in enumerate(encabezados):
            ancho_contenido = max(
                [fuente_encabezado.measure(str(encabezado))]
                + [
                    fuente_celda.measure(str(fila[indice]))
                    for fila in filas_datos
                ]
            )
            anchos.append(min(340, max(110, ancho_contenido + 24)))
        ancho_tabla = sum(anchos) + 40

        def actualizar_region(evento=None):
            region = canvas_tabla.bbox("all")
            if region:
                canvas_tabla.configure(scrollregion=region)

        def ajustar_ancho_canvas(evento):
            canvas_tabla.itemconfigure(
                ventana_tabla, width=max(ancho_tabla, evento.width)
            )

        ventana_tabla = canvas_tabla.create_window(
            (0, 0), window=scrollable_frame, anchor="nw", width=ancho_tabla
        )
        scrollable_frame.bind("<Configure>", actualizar_region)
        canvas_tabla.bind("<Configure>", ajustar_ancho_canvas)
        canvas_tabla.configure(
            yscrollcommand=scrollbar_vertical.set,
            xscrollcommand=scrollbar_horizontal.set,
        )

        scrollable_frame.grid_columnconfigure(0, minsize=anchos[0])
        for indice, encabezado in enumerate(encabezados):
            scrollable_frame.grid_columnconfigure(indice, minsize=anchos[indice])
            tk.Label(
                scrollable_frame, text=encabezado, font=fuente_encabezado,
                fg=self.colores["acento"], bg="#1e293b", anchor="w",
                justify="left", wraplength=max(anchos[indice] - 20, 50),
            ).grid(row=0, column=indice, sticky="nsew", padx=1, pady=1, ipady=9)

        for idx, fila in enumerate(filas_datos):
            bg_color = self.colores["bg_tarjeta"] if idx % 2 == 0 else self.colores["linea_tabla"]
            for indice, val in enumerate(fila):
                texto = str(val)
                try:
                    float(texto.replace(",", ""))
                    es_numero = True
                except ValueError:
                    es_numero = False
                tk.Label(
                    scrollable_frame, text=texto,
                    font=fuente_numero if es_numero else fuente_celda,
                    fg=self.colores["texto_oscuro"], bg=bg_color,
                    anchor="e" if es_numero else "w",
                    justify="right" if es_numero else "left",
                    wraplength=max(anchos[indice] - 20, 50),
                ).grid(
                    row=idx + 1, column=indice, sticky="nsew",
                    padx=10, pady=1, ipady=10,
                )

        scrollbar_vertical.pack(side="right", fill="y", padx=(5, 0))
        scrollbar_horizontal.pack(side="bottom", fill="x")
        canvas_tabla.pack(side="left", fill="both", expand=True)

    def mostrar_resultados(self):
        if not self.resultados:
            self.mostrar_captura()
            return
        self.limpiar_pantalla()
        self._set_active("results")
        self.crear_encabezado(
            "Resultados de las matrículas",
            "Cada fila explica un dato calculado; desplázate horizontalmente si no ves todas las columnas.",
        )

        self.crear_tarjetas_kpi((
            ("Equipo", "4 integrantes", "Matrículas ordenadas para asignar M1–M4.", self.colores["acento"]),
            ("Modelo F", self.modelo["documento_f"]["estado"].capitalize(), "Ecuación térmica original, sin ajustes.", "#f87171" if self.modelo["documento_f"]["estado"] != "completo" else "#34d399"),
            ("Índice G", "No disponible" if self.modelo["documento_g"]["estado"] != "completo" else f"{self.modelo['documento_g']['S']:.5g}", "Solo se calcula con una θ(1) finita.", "#fbbf24"),
        ))
        contenedor = tk.Frame(
            self.content, bg=self.colores["bg_tarjeta"],
            highlightthickness=1, highlightbackground=self.colores["borde"],
        )
        contenedor.pack(fill="both", expand=True, padx=34, pady=(0, 30))

        M, N, P, S = self.resultados
        encabezados = (
            "Dato calculado", "Integrante 1", "Integrante 2", "Integrante 3",
            "Integrante 4", "Qué significa",
        )
        filas = [
            (
                "Matrícula",
                *(
                    "".join(map(str, self.matriculas[f"M{indice}"]))
                    for indice in range(1, 5)
                ),
                "Identificación de cada integrante.",
            ),
            (
                "Número formado (M)",
                *(str(M[f"M{indice}"]) for indice in range(1, 5)),
                SIGNIFICADOS_BASE["M"],
            ),
            (
                "Escala (N)",
                *(f"{N[f'M{indice}']:.7f}" for indice in range(1, 5)),
                SIGNIFICADOS_BASE["N"],
            ),
            (
                "Promedio de dígitos (P)",
                *(f"{P[f'M{indice}']:.7f}" for indice in range(1, 5)),
                SIGNIFICADOS_BASE["P"],
            ),
            (
                "Suma ponderada de dígitos (S base)",
                *(str(S[f"M{indice}"]) for indice in range(1, 5)),
                SIGNIFICADOS_BASE["S"],
            ),
        ]

        self.crear_tabla_resultados(contenedor, filas, encabezados)
        self.agregar_navegacion()

    def mostrar_documento_a(self):
        if not self.resultados:
            self.mostrar_captura()
            return
        self.limpiar_pantalla()
        self._set_active("doc_a")
        self.crear_encabezado(
            "Documento A: Estado operativo",
            "Cada resultado incluye una explicación breve y el valor obtenido.",
        )

        contenedor = tk.Frame(self.content, bg=self.colores["bg_tarjeta"], highlightthickness=1, highlightbackground=self.colores["borde"])
        contenedor.pack(fill="both", expand=True, padx=34, pady=(0, 30))

        vars_a = self.calcular_documento_a()
        filas = crear_filas_descritas(
            (clave, dato["valor"], SIGNIFICADOS_A[clave])
            for clave, dato in vars_a.items()
        )
        self.crear_tabla_resultados(
            contenedor, filas, ("Clave", "Resultado", "Qué significa")
        )
        self.agregar_navegacion()

    def mostrar_documento_b(self):
        self.limpiar_pantalla()
        self._set_active("doc_b")
        self.crear_encabezado(
            "Documento B: Respuesta bajo carga",
            "Los coeficientes describen distintos aportes al modelo de respuesta.",
        )

        contenedor = tk.Frame(self.content, bg=self.colores["bg_tarjeta"], highlightthickness=1, highlightbackground=self.colores["borde"])
        contenedor.pack(fill="both", expand=True, padx=34, pady=(0, 30))

        coefs_b = self.modelo["documento_b"] if self.modelo else calcular_coeficientes_documento_b()
        filas = crear_filas_descritas(
            (clave, coefs_b[clave], SIGNIFICADOS_B[clave])
            for clave in ("Ar", "Br", "Cr", "Dr")
        )
        self.crear_tabla_resultados(
            contenedor, filas, ("Clave", "Resultado", "Qué significa")
        )
        self.agregar_navegacion()

    def mostrar_documento_c(self):
        if not self.resultados:
            self.mostrar_captura()
            return
        self.limpiar_pantalla()
        self._set_active("doc_c")
        self.crear_encabezado(
            "Documento C: Saturación",
            "Valores usados para describir la saturación del sistema, con una explicación de cada uno.",
        )

        contenedor = tk.Frame(self.content, bg=self.colores["bg_tarjeta"], highlightthickness=1, highlightbackground=self.colores["borde"])
        contenedor.pack(fill="both", expand=True, padx=34, pady=(0, 30))

        val_c = self.modelo["documento_c"]
        filas = crear_filas_descritas(
            (clave, val_c[clave], SIGNIFICADOS_C[clave])
            for clave in ("a_s", "b_s", "c_s", "d_s")
        )
        self.crear_tabla_resultados(
            contenedor, filas, ("Clave", "Resultado", "Qué significa")
        )
        self.agregar_navegacion()

    def mostrar_documento_d(self):
        self.limpiar_pantalla()
        self._set_active("doc_d")
        self.crear_encabezado(
            "Documento D: Cambio en el tiempo",
            "Estimaciones obtenidas con tendencias, áreas, búsqueda de raíces y evolución temporal.",
        )
        contenedor = tk.Frame(self.content, bg=self.colores["bg_tarjeta"], highlightthickness=1, highlightbackground=self.colores["borde"])
        contenedor.pack(fill="both", expand=True, padx=34, pady=(0, 30))

        datos_d = self.modelo["documento_d"] if self.modelo else calcular_documento_d_completo()
        filas = crear_filas_descritas((
            ("Pendiente k", datos_d["k"], SIGNIFICADOS_D["k"]),
            ("Valor inicial de la recta a0", datos_d["a0"], SIGNIFICADOS_D["a0"]),
            ("Área estimada m", datos_d["m"], SIGNIFICADOS_D["m"]),
            (
                "Solución n",
                f"{datos_d['n']:.7f} ({datos_d['iteraciones_n']} intentos)",
                SIGNIFICADOS_D["n"],
            ),
            ("Valor estimado p", datos_d["p"], SIGNIFICADOS_D["p"]),
            ("Valor estimado q", datos_d["q"], SIGNIFICADOS_D["q"]),
            ("Balance evaluado L(3)", datos_d["L3"], SIGNIFICADOS_D["L(3)"]),
        ))
        self.crear_tabla_resultados(
            contenedor, filas, ("Resultado", "Valor", "Qué significa")
        )
        self.agregar_navegacion()

    def mostrar_documento_e(self):
        if not self.resultados: return self.mostrar_captura()
        self.limpiar_pantalla()
        self._set_active("doc_e")
        self.crear_encabezado(
            "Documento E: Respuesta del sistema",
            "Estimaciones de cómo cambia la respuesta conforme pasa el tiempo.",
        )
        contenedor = tk.Frame(self.content, bg=self.colores["bg_tarjeta"], highlightthickness=1, highlightbackground=self.colores["borde"])
        contenedor.pack(fill="both", expand=True, padx=34, pady=(0, 30))

        datos_e = self.modelo["documento_e"]

        filas = crear_filas_descritas((
            ("Raíz elegida r_sigma", datos_e["r_sigma"], SIGNIFICADOS_E["r_sigma"]),
            ("Factor de cambio sigma", datos_e["sigma"], SIGNIFICADOS_E["sigma"]),
            ("Nivel objetivo tau", datos_e["tau"], SIGNIFICADOS_E["tau"]),
            ("Aporte total kappa", datos_e["kappa"], SIGNIFICADOS_E["kappa"]),
            (
                "Aporte de referencia omega", datos_e["omega"],
                SIGNIFICADOS_E["omega"],
            ),
            ("Respuesta inicial Y(0)", datos_e["Y0"], SIGNIFICADOS_E["Y(0)"]),
            ("Respuesta final Y(3)", datos_e["Y_final"], SIGNIFICADOS_E["Y(3)"]),
        ))
        self.crear_tabla_resultados(
            contenedor, filas, ("Resultado", "Valor", "Qué significa")
        )
        self.agregar_navegacion()

    def mostrar_documento_f(self):
        if not self.resultados: return self.mostrar_captura()
        self.limpiar_pantalla()
        self._set_active("doc_f")
        self.crear_encabezado(
            "Documento F: Temperatura",
            "Se conserva la ecuación original. La trayectoria puede detenerse si deja de ser finita.",
        )
        datos_f = self.modelo["documento_f"]
        theta0_texto = f"{datos_f['theta0']:.4g} °C"
        self.crear_tarjetas_kpi((
            ("Estado de F", datos_f["estado"].capitalize(), "Sin modificar la ecuación del expediente.", "#f87171" if datos_f["estado"] != "completo" else "#34d399"),
            ("Temperatura inicial", theta0_texto, "Condición propuesta por el equipo.", self.colores["acento"]),
            ("Capacidad logística", f"{datos_f['capacidad_maxima']:.4g}", f"Aporte mínimo Y/ξ: {datos_f['fuente_minima']:.4g}.", "#fbbf24"),
        ))
        contenedor = tk.Frame(self.content, bg=self.colores["bg_tarjeta"], highlightthickness=1, highlightbackground=self.colores["borde"])
        contenedor.pack(fill="both", expand=True, padx=34, pady=(0, 30))

        filas = crear_filas_descritas((
            ("Coeficiente de cambio chi", datos_f["chi"], SIGNIFICADOS_F["chi"]),
            ("Límite elegido psi", datos_f["psi"], SIGNIFICADOS_F["psi"]),
            ("Combinación de valores xi", datos_f["xi"], SIGNIFICADOS_F["xi"]),
            ("Estado de la integración", datos_f["estado"], SIGNIFICADOS_F["estado"]),
            ("Menor aporte Y/xi", datos_f["fuente_minima"], SIGNIFICADOS_F["fuente_minima"]),
            (
                "Mayor aporte original Y/xi", datos_f["fuente_maxima"],
                SIGNIFICADOS_F["Y/xi"],
            ),
            (
                "Capacidad máxima del modelo", datos_f["capacidad_maxima"],
                SIGNIFICADOS_F["capacidad"],
            ),
            (
                "Temperatura inicial theta(0)", datos_f["theta0"],
                SIGNIFICADOS_F["theta(0)"],
            ),
            (
                "Temperatura final theta(3)", datos_f["theta_final"],
                SIGNIFICADOS_F["theta(3)"],
            ),
            (
                "Tiempo donde se detuvo", datos_f["tiempo_fallo"],
                "Tiempo aproximado donde RK4 dejó de avanzar; no aplica si terminó.",
            ),
            (
                "Diagnóstico", datos_f["causa"] or "La integración llegó a t=3.",
                "Explicación del estado de la trayectoria.",
            ),
        ))
        tk.Label(
            contenedor,
            text=(
                f"{datos_f['justificacion_theta0']} "
                "G queda no disponible si no se obtiene una trayectoria finita."
            ),
            bg=self.colores["bg_tarjeta"], fg=self.colores["texto_secundario"],
            font=("Segoe UI", 10), wraplength=950, justify="left",
        ).pack(side="bottom", fill="x", padx=18, pady=8)
        self.crear_tabla_resultados(
            contenedor, filas, ("Resultado", "Valor", "Qué significa")
        )
        self.agregar_navegacion()

    def mostrar_documento_g(self):
        if not self.resultados: return self.mostrar_captura()
        self.limpiar_pantalla()
        self._set_active("doc_g")
        self.crear_encabezado(
            "Documento G: Índice general",
            "Resumen que combina los resultados de respuesta, temperatura y saturación.",
        )
        estado_g = self.modelo["documento_g"]["estado"]
        indice_g = self.modelo["documento_g"].get("S")
        theta0 = self.modelo["documento_f"]["theta0"]
        theta0_texto = f"{theta0:.4g} °C"
        self.crear_tarjetas_kpi((
            ("Estado de G", "No disponible" if indice_g is None else "Calculado", "Depende de una θ(1) finita.", "#fbbf24" if indice_g is None else "#34d399"),
            ("Entrada térmica", self.modelo["documento_f"]["estado"].capitalize(), "Resultado del Documento F original.", "#f87171" if estado_g != "completo" else "#34d399"),
            ("Condición inicial", theta0_texto, "Temperatura inicial propuesta.", self.colores["acento"]),
        ))
        contenedor = tk.Frame(self.content, bg=self.colores["bg_tarjeta"], highlightthickness=1, highlightbackground=self.colores["borde"])
        contenedor.pack(fill="both", expand=True, padx=34, pady=(0, 30))

        datos_g = self.modelo["documento_g"]

        if datos_g["estado"] == "completo":
            datos = (
                ("Respuesta acumulada", datos_g["integral_Y"], SIGNIFICADOS_G["Integral de Y"]),
                ("Temperatura después de 1 unidad", datos_g["theta_1"], SIGNIFICADOS_G["theta(1)"]),
                ("Nivel operativo elegido", datos_g["U_star"], SIGNIFICADOS_G["U*"]),
                ("Índice general (S de G)", datos_g["S"], SIGNIFICADOS_G["S"]),
            )
        else:
            datos = (
                ("Estado", "No disponible", "No se calculó G porque F no proporciona una trayectoria finita para evaluar theta(1)."),
                ("Motivo", datos_g["causa"], "El modelo térmico original se conserva sin ajustes."),
            )
        filas = crear_filas_descritas(datos)
        self.crear_tabla_resultados(
            contenedor, filas, ("Resultado", "Valor", "Qué significa")
        )
        self.agregar_navegacion()

    def dibujar_grafica(
        self, canvas, historial, titulo, color, eje_y, estado=None, tiempo_fallo=None
    ):
        def dibujar(evento=None):
            canvas.delete("all")
            ancho = max(evento.width if evento else canvas.winfo_width(), 200)
            alto = max(evento.height if evento else canvas.winfo_height(), 150)
            izquierda, derecha, arriba, abajo = 64, 18, 38, 48
            x0, x1 = izquierda, ancho - derecha
            y0, y1 = arriba, alto - abajo
            tiempos = [punto[0] for punto in historial]
            valores = [punto[1] for punto in historial]
            t_min, t_max = min(tiempos), max(tiempos)
            etiqueta_eje_y = eje_y
            valores_graficados, recorte_escala = limitar_trayectoria_para_grafica(
                valores, estado
            )
            y_min, y_max = min(valores_graficados), max(valores_graficados)
            if y_min == y_max:
                y_max = y_min + 1
            puntos = []
            for tiempo, valor in zip(tiempos, valores_graficados):
                x = x0 + (tiempo - t_min) * (x1 - x0) / max(t_max - t_min, 1e-12)
                y = y1 - (valor - y_min) * (y1 - y0) / (y_max - y_min)
                puntos.append((x, y))
            canvas.create_text(
                ancho / 2, 17, text=titulo, fill=self.colores["texto_oscuro"],
                font=("Segoe UI", 10, "bold"),
            )
            for indice in range(5):
                proporcion = indice / 4
                x = x0 + proporcion * (x1 - x0)
                y = y1 - proporcion * (y1 - y0)
                tiempo = t_min + proporcion * (t_max - t_min)
                valor = y_min + proporcion * (y_max - y_min)
                canvas.create_line(
                    x, y0, x, y1, fill=self.colores["grid"],
                    dash=(2, 4), width=1,
                )
                canvas.create_line(
                    x0, y, x1, y, fill=self.colores["grid"],
                    dash=(2, 4), width=1,
                )
                canvas.create_text(
                    x, y1 + 16, text=f"{tiempo:.3g}", fill=self.colores["texto_secundario"],
                    font=("Consolas", 8),
                )
                canvas.create_text(
                    x0 - 8, y, text=f"{valor:.3g}", anchor="e",
                    fill=self.colores["texto_secundario"], font=("Consolas", 8),
                )
            canvas.create_line(x0, y1, x1, y1, fill=self.colores["borde"], width=1)
            canvas.create_line(x0, y0, x0, y1, fill=self.colores["borde"], width=1)
            canvas.create_polygon(
                [x0, y1] + [coordenada for punto in puntos for coordenada in punto]
                + [x1, y1],
                fill=color, stipple="gray50", outline="",
            )
            canvas.create_line(
                *[coordenada for punto in puntos for coordenada in punto],
                fill=color, width=2, smooth=False,
            )
            salto = max(1, len(puntos) // 24)
            for indice in range(0, len(puntos), salto):
                x, y = puntos[indice]
                canvas.create_oval(x-2, y-2, x+2, y+2, fill=color, outline="")
            if puntos and (len(puntos) - 1) % salto:
                x, y = puntos[-1]
                canvas.create_oval(x-2, y-2, x+2, y+2, fill=color, outline="")
            canvas.create_text(
                (x0 + x1) / 2, alto - 9, text="Tiempo t",
                fill=self.colores["texto_secundario"], font=("Segoe UI", 8),
            )
            canvas.create_text(
                14, (y0 + y1) / 2, text=etiqueta_eje_y, angle=90,
                fill=self.colores["texto_secundario"], font=("Segoe UI", 8),
            )
            if estado == "divergente" and tiempo_fallo is not None and puntos:
                x, y = puntos[-1]
                canvas.create_line(
                    x, y0, x, y1, fill="#f87171", dash=(3, 3), width=1,
                )
                canvas.create_oval(
                    x-4, y-4, x+4, y+4, fill="#f87171",
                    outline=self.colores["bg_tarjeta"], width=1,
                )
                canvas.create_text(
                    x1-5, y0+9,
                    text=f"{'Escala recortada' if recorte_escala else 'F detenido'} · t ≈ {tiempo_fallo:.3g}",
                    anchor="ne", fill="#fecaca", font=("Segoe UI", 8, "bold"),
                )

        canvas.bind("<Configure>", dibujar)
        self.after_idle(dibujar)

    def mostrar_analisis(self):
        if not self.modelo:
            self.mostrar_captura()
            return
        self.limpiar_pantalla()
        self._set_active("analysis")
        self.crear_encabezado(
            "Análisis de sensibilidad y trayectorias",
            "Compara cuánto cambian las estimaciones al hacer pasos de tiempo más pequeños.",
        )
        estado_f = self.modelo["documento_f"]["estado"]
        indice_g = self.modelo["documento_g"].get("S")
        self.crear_tarjetas_kpi((
            ("Integrador", "Runge-Kutta 4", "Comparación de h = 0.1, 0.05 y 0.025.", self.colores["acento"]),
            ("Modelo F", estado_f.capitalize(), "La trayectoria se conserva hasta el último paso finito.", "#f87171" if estado_f != "completo" else "#34d399"),
            ("Índice G", "No disponible" if indice_g is None else f"{indice_g:.5g}", "Se omite si θ(1) no es finita.", "#fbbf24" if indice_g is None else "#34d399"),
        ))
        panel = tk.Frame(
            self.content, bg=self.colores["bg_tarjeta"],
            highlightthickness=1, highlightbackground=self.colores["borde"],
            padx=14, pady=12,
        )
        panel.pack(fill="both", expand=True, padx=34, pady=(0, 30))
        tabla = tk.Frame(panel, bg=self.colores["bg_tarjeta"])
        tabla.pack(fill="x", padx=12, pady=(8, 2))
        filas = [
            (
                f"{dato['h']:.2f}", f"{dato['Y_final']:.6f}",
                "No disponible" if dato["theta_final"] is None else f"{dato['theta_final']:.6f}",
                "No disponible" if dato["S"] is None else f"{dato['S']:.7f}",
                "—" if dato["diferencia_S"] is None else f"{dato['diferencia_S']:.3g}",
                dato["estado_F"],
                "—" if dato["tiempo_fallo_F"] is None else f"{dato['tiempo_fallo_F']:.4f}",
            )
            for dato in self.modelo["sensibilidad"]["pasos"]
        ]
        self.crear_tabla_resultados(
            tabla, filas, (
                "Paso h", "Respuesta Y(3)", "Temperatura θ(3)",
                "Índice G", "Cambio de G", "Estado F", "Tiempo de fallo",
            )
        )
        detalles_tolerancias = []
        for dato in self.modelo["sensibilidad"]["tolerancias"]:
            indice = (
                "no disponible"
                if dato["S"] is None
                else f"{dato['S']:.6f}"
            )
            detalles_tolerancias.append(
                f"margen {dato['tolerancia']:.0e}: solución {dato['C']:.7f}, "
                f"{dato['iteraciones']} intentos, error restante "
                f"{abs(dato['residuo']):.2g}, índice general {indice}"
            )
        texto_tolerancias = (
            "El tamaño del paso es el intervalo entre cálculos: reducirlo permite "
            "observar si los resultados se mantienen parecidos. "
            "La diferencia del índice compara cada fila con el paso siguiente más pequeño."
            "\n\nSi F diverge, su temperatura final y G se muestran como no disponibles. "
            "\n\nTambién se probó la búsqueda de C con distintos márgenes de error: "
            + " | ".join(detalles_tolerancias)
        )
        tk.Label(
            panel, text=texto_tolerancias, bg=self.colores["bg_tarjeta"],
            fg=self.colores["texto_secundario"], font=("Segoe UI", 9),
            wraplength=950, justify="left",
        ).pack(fill="x", padx=18, pady=(2, 4))
        graficas = tk.Frame(panel, bg=self.colores["bg_tarjeta"])
        graficas.pack(fill="both", expand=True, padx=12, pady=8)
        grafica_y = tk.Canvas(
            graficas, bg=self.colores["bg_tarjeta"], highlightthickness=1,
            highlightbackground=self.colores["borde"], height=265,
        )
        grafica_theta = tk.Canvas(
            graficas, bg=self.colores["bg_tarjeta"], highlightthickness=1,
            highlightbackground=self.colores["borde"], height=265,
        )
        grafica_y.pack(side="left", fill="both", expand=True, padx=4)
        grafica_theta.pack(side="left", fill="both", expand=True, padx=4)
        self.dibujar_grafica(
            grafica_y, self.modelo["documento_e"]["historial"],
            "Documento E · respuesta", self.colores["acento"], "Respuesta Y(t)",
        )
        self.dibujar_grafica(
            grafica_theta, self.modelo["documento_f"]["historial"],
            "Documento F · temperatura", "#f87171", "Temperatura θ(t)",
            self.modelo["documento_f"]["estado"],
            self.modelo["documento_f"]["tiempo_fallo"],
        )
        tk.Label(
            panel,
            text="Los puntos muestran pasos RK4. En F, el extremo se recorta solo para visualizar la trayectoria previa; no se alteran los valores del modelo.",
            bg=self.colores["bg_tarjeta"], fg=self.colores["texto_secundario"],
            font=("Segoe UI", 9), wraplength=950, justify="left",
        ).pack(fill="x", padx=18, pady=(0, 8))
        self.agregar_navegacion()

    def exportar_reporte(self):
        if not self.modelo:
            messagebox.showerror("Sin resultados", "Primero captura cuatro matrículas y calcula el modelo.")
            return
        destino = filedialog.asksaveasfilename(
            title="Guardar informe reproducible",
            defaultextension=".html",
            initialfile="Informe_PIA.html",
            filetypes=(("Informe HTML", "*.html"),),
        )
        if not destino:
            return
        try:
            Path(destino).write_text(
                crear_reporte_html(self.modelo, self.matriculas),
                encoding="utf-8"
            )
        except OSError as error:
            messagebox.showerror("Error al guardar el informe", str(error))
            return
        messagebox.showinfo(
            "Informe guardado",
            "El informe HTML incluye resultados, gráficas y sensibilidad. "
            "Ábrelo en el navegador y usa Imprimir > Guardar como PDF."
        )

    def mostrar_operaciones(self):
        if not self.modelo:
            self.mostrar_captura()
            return
        self.limpiar_pantalla()
        self._set_active("operations")
        self.crear_encabezado(
            "Operaciones e iteraciones detalladas",
            "Resultados intermedios, criterios de paro, supuestos y comparación de pasos."
        )
        contenedor = tk.Frame(
            self.content, bg=self.colores["bg_tarjeta"],
            highlightthickness=1, highlightbackground=self.colores["borde"],
        )
        contenedor.pack(fill="both", expand=True, padx=34, pady=(0, 30))
        texto_area = tk.Text(
            contenedor, bg=self.colores["bg_tarjeta"], fg=self.colores["texto_oscuro"],
            font=("Consolas", 10), relief="flat", padx=20, pady=20, wrap="word",
            insertbackground="white",
        )
        scrollbar = tk.Scrollbar(contenedor, command=texto_area.yview)
        texto_area.configure(yscrollcommand=scrollbar.set)
        texto_area.bind(
            "<MouseWheel>",
            lambda event: (
                texto_area.yview_scroll(-1 if event.delta > 0 else 1, "units"),
                "break",
            )[1],
        )
        texto_area.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

        M, N, P, S_base = self.modelo["resultados"]
        variables_a = self.modelo["variables_a"]
        val_c = self.modelo["documento_c"]
        datos_d = self.modelo["documento_d"]
        val_e = self.modelo["documento_e"]
        val_f = self.modelo["documento_f"]
        val_g = self.modelo["documento_g"]
        coefs_b = self.modelo["documento_b"]
        lineas = [
            "GUÍA RÁPIDA",
            "M es la matrícula escrita como número; N es ese número en una escala pequeña.",
            "P es el promedio de los siete dígitos; S suma los dígitos dando más peso a los últimos.",
            "Respuesta y temperatura son estimaciones del modelo, no mediciones reales.",
            "El índice general S no se clasifica como estable o inestable porque el expediente no da umbrales.",
            "",
            "PARÁMETROS BASE (matrículas ordenadas)",
        ]
        for clave in ("M1", "M2", "M3", "M4"):
            matricula = "".join(map(str, self.matriculas[clave]))
            lineas.append(
                f"{clave}={matricula}: M={M[clave]}, N={N[clave]:.7f}, "
                f"P={P[clave]:.7f}, S={S_base[clave]}"
            )
        lineas.extend(("", "DOCUMENTO A — variables y Newton-Raphson"))
        for clave, dato in variables_a.items():
            lineas.append(f"{clave} [{dato['formula']}] = {dato['valor']:.10g}")
        lineas.append(
            f"C converge cuando |x_nueva-x|<1e-10 y |f(x_nueva)|<1e-10, "
            f"máximo 100 iteraciones; "
            f"convergió en {variables_a['C']['iteraciones']}."
        )
        for paso in variables_a["C"]["historial"]:
            lineas.append(
                f"  iter {paso['iter']}: x={paso['x']:.10g}, "
                f"f(x)={paso['fx']:.4g}, f'(x)={paso['dfx']:.4g}, "
                f"x_nueva={paso['x_nueva']:.10g}"
            )
        lineas.extend(("", "DOCUMENTO B — coeficientes del modelo R(n)"))
        lineas.extend(
            f"{clave}={coefs_b[clave]:.10g}" for clave in ("Ar", "Br", "Cr", "Dr")
        )
        lineas.append(
            "Ar se obtiene aislando las raíces reales de la cúbica y seleccionando "
            "la menor raíz positiva."
        )
        lineas.append(
            "Cr usa mínimos cuadrados sobre (1,4.2),(2,6.1),(4,9.8),(6,13.9); "
            "Dr interpola (1,2.8),(2,4.6),(3,7.9),(4,12.5) en 2.5."
        )
        lineas.extend(("", "DOCUMENTO C — saturación"))
        lineas.extend(
            f"{clave}={val_c[clave]:.10g}" for clave in ("a_s", "b_s", "c_s", "d_s")
        )
        desglose_cs = val_c["desglose_c_s"]
        lineas.append(
            "c_s: Simpson 1/3 compuesto sobre la tabla l(t); "
            f"n={desglose_cs['num_intervalos']} intervalos de paso "
            f"h={desglose_cs['h']:g}, pesos 1, 4, 2, 4, 1."
        )
        for tiempo, valor in zip(desglose_cs["nodos"], desglose_cs["valores"]):
            lineas.append(f"  t={tiempo:g}, l(t)={valor:.4g}")
        lineas.append(
            f"  Σ extremos (l0+l4)={desglose_cs['suma_extremos']:.10g}; "
            f"Σ impares (l1+l3)={desglose_cs['suma_impares']:.10g}; "
            f"Σ pares internos (l2)={desglose_cs['suma_pares']:.10g}"
        )
        lineas.append(
            f"  Σ ponderada = {desglose_cs['suma_extremos']:.10g} "
            f"+ 4·{desglose_cs['suma_impares']:.10g} "
            f"+ 2·{desglose_cs['suma_pares']:.10g} "
            f"= {desglose_cs['suma_ponderada']:.10g}"
        )
        lineas.append(
            f"  c_s = (h/3)·Σ ponderada = ({desglose_cs['h']:g}/3)·"
            f"{desglose_cs['suma_ponderada']:.10g} = {desglose_cs['resultado']:.10g}"
        )
        lineas.append(f"Solución de Gauss-Jordan: x={val_c['x']:.8g}, y={val_c['y']:.8g}, z={val_c['z']:.8g}")
        lineas.append(
            f"Tm (Lagrange en 410)={val_c['Tm']:.8g}; "
            f"Tc (cuadrática evaluada en a)={val_c['Tc']:.8g}"
        )
        lineas.append("Matrices intermedias de Gauss-Jordan:")
        lineas.append("Gauss-Jordan usa pivoteo parcial y tolerancia de pivote 1e-12.")
        for paso in val_c["pasos_gj"]:
            lineas.append(f"  {paso['etapa']}")
            lineas.extend(f"    {fila}" for fila in paso["matriz"])
        lineas.append("Reducción de la matriz para Tc:")
        lineas.extend(
            f"  {paso['etapa']}: {paso['matriz']}"
            for paso in val_c["pasos_gj_tc"]
        )
        lineas.extend(("", "DOCUMENTO D — mínimos cuadrados e integración"))
        lineas.append(
            f"Datos de regresión: {list(zip(datos_d['r_vals'], datos_d['y_vals']))}; "
            f"k={datos_d['k']:.8g}, a0={datos_d['a0']:.8g}"
        )
        lineas.append(
            f"Sumatorias: Σr={datos_d['sum_r']:.8g}, Σy={datos_d['sum_y']:.8g}, "
            f"Σr²={datos_d['sum_r2']:.8g}, Σry={datos_d['sum_ry']:.8g}"
        )
        lineas.extend(
            f"{clave}={datos_d[clave]:.10g}" for clave in ("m", "n", "p", "q")
        )
        lineas.append(
            f"Bisección para n: tolerancia=1e-10, iteraciones={datos_d['iteraciones_n']}"
        )
        for paso in datos_d["historial_n"]:
            lineas.append(
                f"  iter {paso['iter']}: intervalo=[{paso['a']:.10g}, "
                f"{paso['b']:.10g}], x={paso['x']:.10g}, f(x)={paso['fx']:.4g}"
            )
        lineas.append("Runge-Kutta 4 para p: y'=-0.6y+4t, y(0)=1, h=0.1:")
        lineas.extend(
            f"  t={tiempo:.1f}, y={valor:.10g}"
            for tiempo, valor in datos_d["historial_p"]
        )
        lineas.extend(("", "DOCUMENTO E — EDO de rendimiento"))
        lineas.append("Runge-Kutta 4, intervalo [0,3], h=0.1; omega se calcula con Simpson 3/8.")
        lineas.append(
            f"r_sigma={val_e['r_sigma']}, sigma={val_e['sigma']}, "
            f"tau={val_e['tau']:.8g}, kappa={val_e['kappa']:.8g}, "
            f"omega={val_e['omega']:.8g}, h={val_e['h']}"
        )
        lineas.append(f"Dependencia D→E: kappa=L(3)={datos_d['L3']:.10g}.")
        lineas.append(f"Y(0)={val_e['Y0']:.8g}; Y(3)={val_e['Y_final']:.8g}")
        lineas.extend(("", "DOCUMENTO F — EDO térmica"))
        lineas.append("Runge-Kutta 4, intervalo [0,3], h=0.1; Y(t) de E se interpola linealmente entre pasos.")
        lineas.append(
            f"chi={val_f['chi']:.8g}, psi={val_f['psi']}, xi=a+d+g={val_f['xi']:.8g}"
        )
        theta_final = (
            f"{val_f['theta_final']:.8g} °C"
            if val_f["theta_final"] is not None else "No disponible (integración incompleta)"
        )
        lineas.append(
            f"Estado={val_f['estado']}; theta(0)={val_f['theta0']} °C; "
            f"theta(3)={theta_final}"
        )
        lineas.append(f"Justificación: {val_f['justificacion_theta0']}")
        lineas.append(
            f"Ecuación original sin normalización: theta'=chi*theta*(psi-theta)-Y/xi; "
            f"Y/xi en [{val_f['fuente_minima']:.8g}, {val_f['fuente_maxima']:.8g}], "
            f"capacidad logística máxima={val_f['capacidad_maxima']:.8g}."
        )
        if val_f["tiempo_fallo"] is not None:
            lineas.append(
                f"RK4 dejó de avanzar cerca de t={val_f['tiempo_fallo']:.8g}: "
                f"{val_f['causa']}"
            )
        lineas.extend(("", "DOCUMENTO G — índice global"))
        if val_g["estado"] == "completo":
            lineas.append(f"Raíces reales de saturación: {val_g['raices_saturacion']}")
            lineas.append(
                f"Se elige U*={val_g['U_star']:.10g}; integral Y(0,2)="
                f"{val_g['integral_Y']:.10g}; theta(1)={val_g['theta_1']:.10g}; S={val_g['S']:.10g}"
            )
        else:
            lineas.append(f"Índice G no disponible: {val_g['causa']}")
        lineas.extend(("", "ANÁLISIS DE SENSIBILIDAD (Runge-Kutta 4)"))
        for dato in self.modelo["sensibilidad"]["pasos"]:
            diferencia = "—" if dato["diferencia_S"] is None else f"{dato['diferencia_S']:.4g}"
            theta = "No disponible" if dato["theta_final"] is None else f"{dato['theta_final']:.8g}"
            indice = "No disponible" if dato["S"] is None else f"{dato['S']:.8g}"
            fallo = "—" if dato["tiempo_fallo_F"] is None else f"{dato['tiempo_fallo_F']:.5g}"
            lineas.append(
                f"h={dato['h']:.3g}: Y(3)={dato['Y_final']:.8g}, "
                f"theta(3)={theta}, G={indice}, estado F={dato['estado_F']}, "
                f"tiempo de fallo={fallo}, |ΔG|={diferencia}"
            )
        lineas.append("Sensibilidad a tolerancia de Newton-Raphson para C:")
        for dato in self.modelo["sensibilidad"]["tolerancias"]:
            indice = "No disponible" if dato["S"] is None else f"{dato['S']:.10g}"
            lineas.append(
                f"  tol={dato['tolerancia']:.0e}: C={dato['C']:.10g}, "
                f"iteraciones={dato['iteraciones']}, residuo={dato['residuo']:.4g}, "
                f"G={indice}, estado F={dato['estado_F']}"
            )
        lineas.extend((
            "",
            "INTERPRETACIÓN",
            "S es un índice compuesto del expediente. Como no se proporcionan umbrales, no se etiqueta por sí solo como estable o inestable.",
            "",
            "SUPUESTOS DECLARADOS",
            "• Se usa r_sigma=1.2, la menor raíz positiva, para el decaimiento más conservador.",
            "• kappa se calcula como L(3) a partir de los coeficientes del Documento D; omega se integra con Simpson 3/8 en la tabla indicada.",
            "• c_s se calcula con Simpson 1/3 compuesto sobre la tabla l(t) (n=4, h=1); el desglose de sumatorias aparece en esta vista y en el informe HTML.",
            "• Se propone theta(0)=25 °C como línea base ambiente/de reposo previa a la carga y se usa en todos los análisis comparativos.",
            "• Se conserva la ecuación F original sin normalización; cuando Y/xi supera la capacidad logística máxima chi*psi²/4, el aporte negativo domina, la trayectoria pierde finitud y G se declara no disponible.",
            "• Para U* se elige la menor raíz real positiva de la ecuación de saturación.",
        ))
        texto_area.tag_configure(
            "section", foreground=self.colores["acento"],
            font=("Segoe UI", 11, "bold"), spacing1=12, spacing3=4,
        )
        texto_area.tag_configure(
            "formula", foreground="#c4b5fd", font=("Consolas", 10),
        )
        texto_area.tag_configure(
            "result", foreground="#a7f3d0", font=("Consolas", 10, "bold"),
        )
        texto_area.tag_configure(
            "warning", foreground="#fca5a5", font=("Segoe UI", 10, "bold"),
            background="#2a1820",
        )
        texto_area.tag_configure(
            "muted", foreground=self.colores["texto_secundario"],
        )
        texto_area.insert("1.0", "\n".join(lineas))
        for indice, linea in enumerate(lineas, start=1):
            rango = f"{indice}.0", f"{indice}.end"
            linea_limpia = linea.strip()
            if linea_limpia.startswith((
                "DOCUMENTO ", "ANÁLISIS DE SENSIBILIDAD",
                "SUPUESTOS DECLARADOS", "INTERPRETACIÓN",
                "GUÍA RÁPIDA", "PARÁMETROS BASE",
            )):
                texto_area.tag_add("section", *rango)
            elif any(clave in linea_limpia.lower() for clave in (
                "no disponible", "divergente", "dejó de avanzar", "no es válido",
            )):
                texto_area.tag_add("warning", *rango)
            elif linea_limpia.startswith((
                "C converge", "r_sigma=", "Y(0)=", "Estado=",
                "Se elige U*=", "Ecuación original",
            )) or " y(" in linea_limpia.lower() or " theta'" in linea_limpia.lower():
                texto_area.tag_add("result", *rango)
            elif linea.startswith("  iter ") or linea.startswith("    "):
                texto_area.tag_add("muted", *rango)
        texto_area.configure(state="disabled")
        self.agregar_navegacion()


if __name__ == "__main__":
    Aplicacion().mainloop()