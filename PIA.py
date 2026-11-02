import tkinter as tk
import tkinter.font as tkfont
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
    # ---------- Ar: raíz positiva menor de w^3 - 2.1w^2 + 0.9w - 0.08 = 0 ----------
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
    Dr = interpolacion_lagrange(2.5, nodos_dr, valores_dr)

    return {
        "Ar": Ar,
        "Br": Br,
        "Cr": Cr,
        "Dr": Dr,
        "historial_ar": historial_ar,
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

# ==========================================
# NUEVOS DOCUMENTOS (D, E, F, G)
# ==========================================
def calcular_documento_d_completo():
    datos_regresion = calcular_documento_d()
    m_val = simpson_1_3(lambda u: 3*u**2 + 1.5*u + 2, 0, 2, 2)
    n_val, historial_n = metodo_biseccion_historial(
        lambda x: x**3 - 7*x + 6, 1 + 1e-8, 3, tolerancia=1e-10
    )
    historial_p = euler_adelante_historial(
        lambda t, y: -0.6*y + 4*t, 0, 1, 1.5, 0.1
    )
    p_val = historial_p[-1][1]
    q_val = interpolacion_lagrange(4.5, (2, 4, 6, 8), (10.5, 18.7, 31.1, 47.3))

    return {
        **datos_regresion, "m": m_val, "n": n_val, "p": p_val, "q": q_val,
        "historial_n": historial_n, "iteraciones_n": len(historial_n),
        "historial_p": historial_p,
    }

def calcular_documento_e_completo(val_c, h=0.1):
    Ar, Tm, x = val_c["Ar"], val_c["Tm"], val_c["x"]

    # Se elige la raíz menor positiva para usar el decaimiento más conservador.
    r_sigma = 1.2
    sigma = r_sigma / 10.0
    tau = Tm - Ar
    kappa = simpson_3_8((1.5, 2.6, 4.8, 8.3), 1.0)

    R_1_8 = Ar*(1.8**3) + val_c["Br"]*(1.8**2) + val_c["Cr"]*1.8 + val_c["Dr"]
    Y0 = R_1_8 + x
    omega = kappa

    def edo(t, Y):
        if t < 0:
            raise ValueError("El tiempo del modelo E debe ser no negativo.")
        return -sigma*(Y - tau) + kappa/omega

    historial = euler_adelante_historial(edo, 0, Y0, 3.0, h)

    return {
        "r_sigma": r_sigma, "sigma": sigma, "tau": tau,
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
    psi = 8.0
    xi = variables_a["a"]["valor"] + variables_a["d"]["valor"] + variables_a["g"]["valor"]
    if abs(xi) < 1e-14:
        raise ValueError("xi = a + d + g no puede ser cero.")

    theta0 = 25.0
    capacidad_maxima = chi * psi**2 / 4
    fuente_maxima = max(
        interpolar_historial(val_e["historial"], t) / xi
        for t, _ in val_e["historial"]
    )
    fraccion_capacidad = 0.8
    escala_fuente = max(1.0, fuente_maxima / (fraccion_capacidad * capacidad_maxima))
    edo = lambda t, T: (
        chi*T*(psi - T)
        - interpolar_historial(val_e["historial"], t)/(xi * escala_fuente)
    )
    historial = euler_adelante_historial(edo, 0, theta0, 3.0, h)

    return {
        "chi": chi, "psi": psi, "xi": xi, "theta0": theta0,
        "theta_final": historial[-1][1], "historial": historial,
        "h": h, "escala_fuente": escala_fuente,
        "fuente_maxima": fuente_maxima,
        "capacidad_maxima": capacidad_maxima,
        "fraccion_capacidad": fraccion_capacidad,
        "justificacion_theta0": (
            "Se propone 25 °C como temperatura ambiente/de reposo típica "
            "antes de iniciar la carga de trabajo."
        ),
    }

def calcular_documento_g_completo(val_e, val_f):
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
        "integral_Y": integral_Y, "theta_1": theta_1,
        "U_star": U_star, "raices_saturacion": raices,
        "S": S,
    }

def calcular_sensibilidad(
    resultados, variables_a, val_c, val_g, pasos=(0.1, 0.05, 0.025)
):
    sensibilidad_pasos = []
    for h in pasos:
        val_e = calcular_documento_e_completo(val_c, h)
        val_f = calcular_documento_f_completo(variables_a, val_e, h)
        val_g = calcular_documento_g_completo(val_e, val_f)
        sensibilidad_pasos.append({
            "h": h, "Y_final": val_e["Y_final"],
            "theta_final": val_f["theta_final"], "S": val_g["S"],
        })
    for indice, resultado in enumerate(sensibilidad_pasos):
        resultado["diferencia_S"] = (
            abs(resultado["S"] - sensibilidad_pasos[indice + 1]["S"])
            if indice + 1 < len(sensibilidad_pasos) else None
        )

    P, S_base = resultados[2], resultados[3]
    sensibilidad_tolerancias = []
    for tolerancia in (1e-6, 1e-8, 1e-10):
        variables_tolerancia = calcular_variables_documento_a(
            resultados, tolerancia_newton=tolerancia
        )
        val_c_tolerancia = calcular_documento_c(variables_tolerancia)
        val_e_tolerancia = calcular_documento_e_completo(val_c_tolerancia)
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
            "S": val_g_tolerancia["S"],
            "diferencia_S": abs(val_g_tolerancia["S"] - val_g["S"]),
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
    val_c = calcular_documento_c(variables_a)
    datos_d = calcular_documento_d_completo()
    val_e = calcular_documento_e_completo(val_c)
    val_f = calcular_documento_f_completo(variables_a, val_e)
    val_g = calcular_documento_g_completo(val_e, val_f)
    sensibilidad = calcular_sensibilidad(
        resultados, variables_a, val_c, val_g
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
    "c_s": "Acumulado de las mediciones entre los tiempos 0 y 4.",
    "d_s": "Combinación de las estimaciones Tm y Tc calculadas con los datos del expediente.",
}

SIGNIFICADOS_D = {
    "k": "Cambio promedio de las mediciones por cada unidad de r.",
    "a0": "Valor de y que predice la recta cuando r vale cero.",
    "m": "Área estimada bajo la curva del modelo entre 0 y 2.",
    "n": "Valor entre 1 y 3 que hace que la ecuación se acerque a cero.",
    "p": "Valor estimado de y cuando el tiempo llega a 1.5.",
    "q": "Valor estimado por la tabla cuando r vale 4.5.",
}

SIGNIFICADOS_E = {
    "r_sigma": "Raíz elegida para definir el comportamiento de la respuesta.",
    "sigma": "Factor que determina qué tan rápido cambia la respuesta.",
    "tau": "Nivel objetivo calculado con los datos recuperados.",
    "kappa": "Aporte total estimado usando cuatro mediciones.",
    "Y(0)": "Respuesta calculada al inicio del intervalo.",
    "Y(3)": "Respuesta estimada después de tres unidades de tiempo.",
}

SIGNIFICADOS_F = {
    "chi": "Qué tan rápido cambia la temperatura según las mediciones disponibles.",
    "psi": "Valor elegido para definir el límite de temperatura del modelo.",
    "xi": "Combinación de tres valores operativos calculados en el Documento A.",
    "lambda": "Factor que reduce el efecto térmico para evitar que el cálculo diverja.",
    "Y/xi": "Intensidad original del efecto de Y antes de aplicar la corrección.",
    "capacidad": "Límite máximo del término logístico de la ecuación térmica.",
    "theta(0)": "Temperatura inicial propuesta: 25 °C antes de iniciar la carga.",
    "theta(3)": "Temperatura estimada al final de las tres unidades de tiempo.",
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
            valor if isinstance(valor, str) else format(valor, formato),
            significado,
        )
        for nombre, valor, significado in datos
    ]

def crear_svg_trayectoria(historial, titulo, color):
    ancho, alto = 640, 250
    margen = 45
    tiempos = [punto[0] for punto in historial]
    valores = [punto[1] for punto in historial]
    t_min, t_max = min(tiempos), max(tiempos)
    y_min, y_max = min(valores), max(valores)
    if y_min == y_max:
        y_max = y_min + 1
    puntos = []
    for t, valor in historial:
        x = margen + (t - t_min) * (ancho - 2*margen) / max(t_max - t_min, 1e-12)
        y = alto - margen - (valor - y_min) * (alto - 2*margen) / (y_max - y_min)
        puntos.append(f"{x:.2f},{y:.2f}")
    return (
        f'<svg viewBox="0 0 {ancho} {alto}" role="img" aria-label="{escape(titulo)}">'
        f'<text x="{ancho/2}" y="22" text-anchor="middle">{escape(titulo)}</text>'
        f'<line x1="{margen}" y1="{alto-margen}" x2="{ancho-margen}" y2="{alto-margen}" stroke="#64748b"/>'
        f'<line x1="{margen}" y1="{margen}" x2="{margen}" y2="{alto-margen}" stroke="#64748b"/>'
        f'<polyline fill="none" stroke="{color}" stroke-width="3" points="{" ".join(puntos)}"/>'
        f'<text x="{margen}" y="{alto-10}">{t_min:.2f}</text>'
        f'<text x="{ancho-margen}" y="{alto-10}" text-anchor="end">{t_max:.2f} t</text>'
        f'<text x="{margen-5}" y="{margen-8}" text-anchor="end">{y_max:.3g}</text>'
        f'<text x="{margen-5}" y="{alto-margen}" text-anchor="end">{y_min:.3g}</text>'
        "</svg>"
    )

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
        ("E", "r_sigma", SIGNIFICADOS_E["r_sigma"], e["r_sigma"]),
        ("E", "sigma", SIGNIFICADOS_E["sigma"], e["sigma"]),
        ("E", "tau", SIGNIFICADOS_E["tau"], e["tau"]),
        ("E", "kappa", SIGNIFICADOS_E["kappa"], e["kappa"]),
        ("E", "Y(0)", SIGNIFICADOS_E["Y(0)"], e["Y0"]),
        ("E", "Y(3)", SIGNIFICADOS_E["Y(3)"], e["Y_final"]),
        ("F", "chi", SIGNIFICADOS_F["chi"], f["chi"]),
        ("F", "psi", SIGNIFICADOS_F["psi"], f["psi"]),
        ("F", "xi", SIGNIFICADOS_F["xi"], f["xi"]),
        (
            "F", "Factor de normalización térmica (λ)",
            SIGNIFICADOS_F["lambda"], f["escala_fuente"],
        ),
        ("F", "Máximo original Y/xi", SIGNIFICADOS_F["Y/xi"], f["fuente_maxima"]),
        ("F", "Capacidad máxima", SIGNIFICADOS_F["capacidad"], f["capacidad_maxima"]),
        ("F", "Temperatura inicial θ(0)", SIGNIFICADOS_F["theta(0)"], f["theta0"]),
        ("F", "Temperatura final θ(3)", SIGNIFICADOS_F["theta(3)"], f["theta_final"]),
        ("G", "Integral de Y", SIGNIFICADOS_G["Integral de Y"], g["integral_Y"]),
        ("G", "Temperatura θ(1)", SIGNIFICADOS_G["theta(1)"], g["theta_1"]),
        ("G", "Nivel operativo U*", SIGNIFICADOS_G["U*"], g["U_star"]),
        ("G", "Índice general S", SIGNIFICADOS_G["S"], g["S"]),
    )
    filas_modelos = "".join(
        f"<tr><td>Documento {escape(documento)}</td><td>{escape(nombre)}</td>"
        f"<td>{escape(descripcion)}</td><td>{valor:.10g}</td></tr>"
        for documento, nombre, descripcion, valor in resultados_documentos
    )
    filas_sensibilidad = []
    for dato in modelo["sensibilidad"]["pasos"]:
        diferencia = "—" if dato["diferencia_S"] is None else f"{dato['diferencia_S']:.3g}"
        filas_sensibilidad.append(
            f"<tr><td>{dato['h']:.3g}</td><td>{dato['Y_final']:.8g}</td>"
            f"<td>{dato['theta_final']:.8g}</td><td>{dato['S']:.8g}</td>"
            f"<td>{diferencia}</td></tr>"
        )
    filas_sensibilidad = "".join(filas_sensibilidad)
    filas_tolerancias = "".join(
        f"<tr><td>{dato['tolerancia']:.0e}</td><td>{dato['C']:.10g}</td>"
        f"<td>{dato['iteraciones']}</td><td>{abs(dato['residuo']):.3g}</td>"
        f"<td>{dato['S']:.10g}</td><td>{dato['diferencia_S']:.3g}</td></tr>"
        for dato in modelo["sensibilidad"]["tolerancias"]
    )
    residuo_n = d["n"]**3 - 7*d["n"] + 6
    U_star = g["U_star"]
    residuo_u = c["a_s"]*U_star**3 - c["b_s"]*U_star**2 + c["c_s"]*U_star - c["d_s"]
    C = variables["C"]["valor"]
    residuo_c = C**3 - (P["M1"] + P["M2"])*C**2 + S_base["M3"]*C - 5
    codigo_fuente = escape(Path(__file__).read_text(encoding="utf-8"))
    return f"""<!doctype html>
<html lang="es"><meta charset="utf-8"><title>Informe técnico PIA</title>
<style>
body{{font:16px/1.5 Arial,sans-serif;max-width:1000px;margin:2rem auto;padding:0 1rem;color:#172033}}
h1,h2{{color:#123b68}} table{{border-collapse:collapse;width:100%;margin:1rem 0 2rem}}
th,td{{border:1px solid #9aa8b8;padding:.45rem;text-align:left}} th{{background:#e8eef5}}
svg{{width:min(100%,700px);border:1px solid #cbd5e1;margin:1rem 0}}
pre{{white-space:pre-wrap;overflow-wrap:anywhere;font-size:8pt;line-height:1.25}}
@media print{{body{{max-width:none}}}}
</style>
<h1>Informe técnico — Proyecto Integrador de Aprendizaje</h1>
<p>Este informe reúne los resultados del modelo en lenguaje sencillo. Las fórmulas y los pasos numéricos se conservan en las secciones técnicas y en el código incluido al final.</p>
<h2>Datos particulares</h2>
<p><b>Número formado (M):</b> matrícula escrita como un número. <b>Escala (N):</b> M dividido entre diez millones. <b>Promedio (P):</b> promedio de sus siete dígitos. <b>Suma ponderada (S base):</b> suma que da más peso a los dígitos ubicados al final. Esta S base es distinta del índice general S del Documento G.</p>
<table><thead><tr><th>Integrante</th><th>Matrícula</th><th>Número formado (M)</th><th>Escala (N)</th><th>Promedio de dígitos (P)</th><th>Suma ponderada (S)</th></tr></thead><tbody>{filas_base}</tbody></table>
<h2>Documento A — Variables operativas</h2>
<p>La columna «En palabras simples» explica para qué sirve cada resultado; la columna técnica conserva la operación utilizada.</p>
<table><thead><tr><th>Clave</th><th>En palabras simples</th><th>Operación técnica</th><th>Resultado</th></tr></thead><tbody>{filas_a}</tbody></table>
<h2>Documentos B–G — Resultados explicados</h2>
<table><thead><tr><th>Documento</th><th>Clave del modelo</th><th>Qué significa</th><th>Resultado</th></tr></thead><tbody>{filas_modelos}</tbody></table>
<h2>Trayectorias numéricas</h2>
<p>En ambas gráficas, el eje horizontal representa el tiempo y el vertical el valor estimado. La línea muestra cómo cambia la respuesta o la temperatura durante el intervalo.</p>
{crear_svg_trayectoria(e["historial"], "Documento E: respuesta del sistema a lo largo del tiempo", "#2563eb")}
{crear_svg_trayectoria(f["historial"], "Documento F: temperatura estimada a lo largo del tiempo", "#dc2626")}
<h2>Sensibilidad al paso de Euler</h2>
<p>Se repitieron los modelos E–G con tres tamaños de paso. Un paso menor divide el tiempo en intervalos más pequeños. La última columna compara el índice general con el resultado del siguiente paso más pequeño; una diferencia menor sugiere que el cálculo cambia menos al refinarlo.</p>
<table><thead><tr><th>Tamaño del paso</th><th>Respuesta al final Y(3)</th><th>Temperatura al final θ(3)</th><th>Índice general S</th><th>Diferencia del índice frente al siguiente paso</th></tr></thead><tbody>{filas_sensibilidad}</tbody></table>
<h3>Sensibilidad a la tolerancia de Newton-Raphson en C</h3>
<p>La tolerancia es el margen de error usado para decidir cuándo detener la búsqueda de C. El residuo indica qué tan cerca queda la solución de cumplir su ecuación; «iteraciones» cuenta los intentos de cálculo.</p>
<table><thead><tr><th>Margen de error permitido</th><th>Solución C</th><th>Intentos de cálculo</th><th>Error restante en la ecuación</th><th>Índice general S</th><th>Diferencia de S frente al cálculo de referencia</th></tr></thead><tbody>{filas_tolerancias}</tbody></table>
<h2>Métodos, criterios y supuestos</h2>
<ul>
<li>Todos los métodos numéricos están implementados en el código: Newton-Raphson, bisección, interpolación de Lagrange, mínimos cuadrados, Gauss-Jordan con pivoteo, trapecio compuesto, Simpson 1/3 y 3/8, y Euler hacia adelante.</li>
<li>Para el modelo E se selecciona r_sigma = 1.2 entre las raíces admisibles 1.2 y 4.0; se usa la raíz menor para adoptar el decaimiento más conservador.</li>
<li>El expediente no proporciona un valor independiente para omega. Se declara el supuesto reproducible omega = kappa = L(3), por lo que kappa/omega = 1; debe confirmarse con el docente.</li>
<li>En F se propone theta(0) = 25 °C como temperatura ambiente/de reposo previa a la carga de trabajo, conforme a la justificación del equipo.</li>
<li>Como el signo de -Y/xi es negativo, ese término reduce theta en la ecuación transcrita; no debe describirse como una generación de calor que eleva la temperatura.</li>
<li>La ecuación F transcrita produce divergencia para entradas típicas porque Y/xi excede la capacidad máxima chi*psi²/4. Para generar una trayectoria finita, el programa usa theta' = chi*theta*(psi-theta) - Y/(xi*lambda), con lambda=max(1, max(Y/xi)/(0.8*chi*psi²/4)). El margen de 20 % es una normalización añadida por el equipo, no un dato del expediente.</li>
<li>Para U* se resuelve la ecuación cúbica de saturación y se elige la menor raíz real positiva. El criterio prioriza el nivel operativo positivo más bajo.</li>
<li>El análisis de sensibilidad compara el índice S entre pasos consecutivos. Los valores son aproximaciones dependientes del paso de Euler.</li>
</ul>
<h2>Interpretación técnica</h2>
<p>El valor S es el índice compuesto definido por el expediente, calculado con las aproximaciones indicadas. El PDF no proporciona umbrales para clasificar el sistema como estable/inestable; por ello, el número debe reportarse junto con sus supuestos y sensibilidad, sin asignarle una categoría física no justificada.</p>
<h2>Validaciones numéricas</h2>
<table><thead><tr><th>Comprobación</th><th>Resultado</th></tr></thead><tbody>
<tr><td>Pendiente k con la tabla del documento D</td><td>{d["k"]:.10g} (referencia manual: 1.71)</td></tr>
<tr><td>Integral m, polinomio integrado exactamente por Simpson</td><td>{d["m"]:.10g} (valor analítico: 15)</td></tr>
<tr><td>Residuo de la raíz n de D</td><td>{residuo_n:.3g}</td></tr>
<tr><td>Residuo de Newton-Raphson para C</td><td>{residuo_c:.3g}</td></tr>
<tr><td>Residuo de la raíz U* de saturación</td><td>{residuo_u:.3g}</td></tr>
</tbody></table>
<h2>Conclusión</h2>
<p>Se reconstruyeron A–G y se hicieron reproducibles desde la captura de las cuatro matrículas. Las decisiones no fijadas en el expediente y la corrección térmica se declaran de manera visible; el índice debe interpretarse junto con su sensibilidad y no como una clasificación absoluta.</p>
<h2>Reproducibilidad y ejecución</h2>
<p>Programa: PIA.py. Requiere Python 3.10 o posterior y Tkinter; no requiere paquetes externos. Ejecutar con <code>python PIA.py</code>, introducir cuatro matrículas de siete dígitos y usar «Exportar informe» para guardar este documento. Desde el navegador se puede imprimir a PDF.</p>
<h2>Código fuente completo</h2><p>El siguiente código es el archivo ejecutable PIA.py usado para generar el informe.</p>
<pre><code>{codigo_fuente}</code></pre>
</html>"""

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
        self.modelo = None
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

        left_container = tk.Frame(barra, bg=self.colores["barra_superior"])
        left_container.pack(side="left", padx=10)
        estilo_boton = {
            "bg": self.colores["barra_superior"], "fg": "#e5e7eb",
            "activebackground": "#1f2937", "activeforeground": "white",
            "relief": "flat", "padx": 9, "pady": 8,
            "font": ("Segoe UI", 10, "bold"), "cursor": "hand2",
        }
        tk.Button(
            left_container, text="Resultados", command=self.mostrar_resultados,
            **estilo_boton
        ).pack(side="left", padx=2, pady=10)
        tk.Button(
            left_container, text="Análisis y gráficas", command=self.mostrar_analisis,
            **estilo_boton
        ).pack(side="left", padx=2, pady=10)

        menu_documentos = tk.Menubutton(left_container, text="Documentos ▾", **estilo_boton)
        submenu = tk.Menu(menu_documentos, tearoff=False)
        for texto, comando in (
            ("Documento A", self.mostrar_documento_a),
            ("Documento B", self.mostrar_documento_b),
            ("Documento C", self.mostrar_documento_c),
            ("Documento D", self.mostrar_documento_d),
            ("Documento E", self.mostrar_documento_e),
            ("Documento F", self.mostrar_documento_f),
            ("Documento G", self.mostrar_documento_g),
            ("Operaciones detalladas", self.mostrar_operaciones),
        ):
            submenu.add_command(label=texto, command=comando)
        menu_documentos.configure(menu=submenu)
        menu_documentos.pack(side="left", padx=2, pady=10)
        tk.Button(
            left_container, text="Exportar informe", command=self.exportar_reporte,
            **estilo_boton
        ).pack(side="left", padx=2, pady=10)

        right_container = tk.Frame(barra, bg=self.colores["barra_superior"])
        right_container.pack(side="right", padx=10)

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
        canvas_tabla = tk.Canvas(contenedor, bg=self.colores["bg_tarjeta"], highlightthickness=0)
        scrollbar_vertical = tk.Scrollbar(contenedor, orient="vertical", command=canvas_tabla.yview)
        scrollbar_horizontal = tk.Scrollbar(
            contenedor, orient="horizontal", command=canvas_tabla.xview
        )
        scrollable_frame = tk.Frame(canvas_tabla, bg=self.colores["bg_tarjeta"])
        fuente_celda = tkfont.Font(font=("Segoe UI", 10))
        fuente_encabezado = tkfont.Font(font=("Segoe UI", 11, "bold"))
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
                fg="#f3f4f6", bg="#1e3a8a", anchor="nw",
                justify="left", wraplength=max(anchos[indice] - 20, 50),
            ).grid(row=0, column=indice, sticky="nsew", padx=10, pady=(10, 5))

        for idx, fila in enumerate(filas_datos):
            bg_color = self.colores["bg_tarjeta"] if idx % 2 == 0 else self.colores["linea_tabla"]
            for indice, val in enumerate(fila):
                tk.Label(
                    scrollable_frame, text=str(val), font=fuente_celda,
                    fg=self.colores["texto_oscuro"], bg=bg_color, anchor="nw",
                    justify="left", wraplength=max(anchos[indice] - 20, 50),
                ).grid(
                    row=idx + 1, column=indice, sticky="nsew",
                    padx=10, pady=12,
                )

        scrollbar_vertical.pack(side="right", fill="y")
        scrollbar_horizontal.pack(side="bottom", fill="x")
        canvas_tabla.pack(side="left", fill="both", expand=True)

    def mostrar_resultados(self):
        if not self.resultados:
            self.mostrar_captura()
            return
        self.limpiar_pantalla()
        self.configure(bg=self.colores["bg_principal"])
        self.crear_encabezado(
            "Resultados de las matrículas",
            "Cada fila explica un dato calculado; desplázate horizontalmente si no ves todas las columnas.",
        )

        contenedor = tk.Frame(self, bg=self.colores["bg_tarjeta"], bd=1, relief="solid")
        contenedor.place(relx=0.05, rely=0.18, relwidth=0.9, relheight=0.68)

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
        self.configure(bg=self.colores["bg_principal"])
        self.crear_encabezado(
            "Documento A: Estado operativo",
            "Cada resultado incluye una explicación breve y el valor obtenido.",
        )

        contenedor = tk.Frame(self, bg=self.colores["bg_tarjeta"], bd=1, relief="solid")
        contenedor.place(relx=0.05, rely=0.18, relwidth=0.9, relheight=0.68)

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
        self.configure(bg=self.colores["bg_principal"])
        self.crear_encabezado(
            "Documento B: Respuesta bajo carga",
            "Los coeficientes describen distintos aportes al modelo de respuesta.",
        )

        contenedor = tk.Frame(self, bg=self.colores["bg_tarjeta"], bd=1, relief="solid")
        contenedor.place(relx=0.05, rely=0.18, relwidth=0.9, relheight=0.68)

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
        self.configure(bg=self.colores["bg_principal"])
        self.crear_encabezado(
            "Documento C: Saturación",
            "Valores usados para describir la saturación del sistema, con una explicación de cada uno.",
        )

        contenedor = tk.Frame(self, bg=self.colores["bg_tarjeta"], bd=1, relief="solid")
        contenedor.place(relx=0.05, rely=0.18, relwidth=0.9, relheight=0.68)

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
        self.configure(bg=self.colores["bg_principal"])
        self.crear_encabezado(
            "Documento D: Cambio en el tiempo",
            "Estimaciones obtenidas con tendencias, áreas, búsqueda de raíces y evolución temporal.",
        )
        contenedor = tk.Frame(self, bg=self.colores["bg_tarjeta"], bd=1, relief="solid")
        contenedor.place(relx=0.05, rely=0.18, relwidth=0.9, relheight=0.68)

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
        ))
        self.crear_tabla_resultados(
            contenedor, filas, ("Resultado", "Valor", "Qué significa")
        )
        self.agregar_navegacion()

    def mostrar_documento_e(self):
        if not self.resultados: return self.mostrar_captura()
        self.limpiar_pantalla()
        self.crear_encabezado(
            "Documento E: Respuesta del sistema",
            "Estimaciones de cómo cambia la respuesta conforme pasa el tiempo.",
        )
        contenedor = tk.Frame(self, bg=self.colores["bg_tarjeta"], bd=1, relief="solid")
        contenedor.place(relx=0.05, rely=0.18, relwidth=0.9, relheight=0.68)

        datos_e = self.modelo["documento_e"]

        filas = crear_filas_descritas((
            ("Raíz elegida r_sigma", datos_e["r_sigma"], SIGNIFICADOS_E["r_sigma"]),
            ("Factor de cambio sigma", datos_e["sigma"], SIGNIFICADOS_E["sigma"]),
            ("Nivel objetivo tau", datos_e["tau"], SIGNIFICADOS_E["tau"]),
            ("Aporte total kappa", datos_e["kappa"], SIGNIFICADOS_E["kappa"]),
            (
                "Aporte de referencia omega", datos_e["omega"],
                "Valor supuesto igual a kappa porque el documento no proporciona otro.",
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
        self.crear_encabezado(
            "Documento F: Temperatura",
            "La temperatura inicial de 25 °C es una propuesta del equipo; el factor lambda es una corrección declarada.",
        )
        contenedor = tk.Frame(self, bg=self.colores["bg_tarjeta"], bd=1, relief="solid")
        contenedor.place(relx=0.05, rely=0.18, relwidth=0.9, relheight=0.68)

        datos_f = self.modelo["documento_f"]

        filas = crear_filas_descritas((
            ("Coeficiente de cambio chi", datos_f["chi"], SIGNIFICADOS_F["chi"]),
            ("Límite elegido psi", datos_f["psi"], SIGNIFICADOS_F["psi"]),
            ("Combinación de valores xi", datos_f["xi"], SIGNIFICADOS_F["xi"]),
            (
                "Factor de corrección lambda", datos_f["escala_fuente"],
                SIGNIFICADOS_F["lambda"],
            ),
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
        ))
        self.crear_tabla_resultados(
            contenedor, filas, ("Resultado", "Valor", "Qué significa")
        )
        tk.Label(
            contenedor,
            text=(
                f"Condición inicial justificada: {datos_f['justificacion_theta0']} "
                "La normalización térmica es una corrección declarada; consulta Operaciones detalladas."
            ),
            bg=self.colores["bg_tarjeta"], fg=self.colores["texto_secundario"],
            font=("Segoe UI", 10), wraplength=950, justify="left",
        ).pack(side="bottom", fill="x", padx=18, pady=8)
        self.agregar_navegacion()

    def mostrar_documento_g(self):
        if not self.resultados: return self.mostrar_captura()
        self.limpiar_pantalla()
        self.crear_encabezado(
            "Documento G: Índice general",
            "Resumen que combina los resultados de respuesta, temperatura y saturación.",
        )
        contenedor = tk.Frame(self, bg=self.colores["bg_tarjeta"], bd=1, relief="solid")
        contenedor.place(relx=0.05, rely=0.18, relwidth=0.9, relheight=0.68)

        datos_g = self.modelo["documento_g"]

        filas = crear_filas_descritas((
            ("Respuesta acumulada", datos_g["integral_Y"], SIGNIFICADOS_G["Integral de Y"]),
            ("Temperatura después de 1 unidad", datos_g["theta_1"], SIGNIFICADOS_G["theta(1)"]),
            ("Nivel operativo elegido", datos_g["U_star"], SIGNIFICADOS_G["U*"]),
            ("Índice general (S de G)", datos_g["S"], SIGNIFICADOS_G["S"]),
        ))
        self.crear_tabla_resultados(
            contenedor, filas, ("Resultado", "Valor", "Qué significa")
        )
        self.agregar_navegacion()

    def dibujar_grafica(self, canvas, historial, titulo, color):
        def dibujar(evento=None):
            canvas.delete("all")
            ancho = max(evento.width if evento else canvas.winfo_width(), 200)
            alto = max(evento.height if evento else canvas.winfo_height(), 150)
            margen = 48
            tiempos = [punto[0] for punto in historial]
            valores = [punto[1] for punto in historial]
            t_min, t_max = min(tiempos), max(tiempos)
            y_min, y_max = min(valores), max(valores)
            if y_min == y_max:
                y_max = y_min + 1
            puntos = []
            for tiempo, valor in historial:
                x = margen + (tiempo - t_min) * (ancho - 2*margen) / max(t_max - t_min, 1e-12)
                y = alto - margen - (valor - y_min) * (alto - 2*margen) / (y_max - y_min)
                puntos.extend((x, y))
            canvas.create_text(ancho/2, 18, text=titulo, fill=self.colores["texto_oscuro"], font=("Segoe UI", 11, "bold"))
            canvas.create_line(margen, alto-margen, ancho-margen, alto-margen, fill=self.colores["texto_secundario"])
            canvas.create_line(margen, margen, margen, alto-margen, fill=self.colores["texto_secundario"])
            canvas.create_line(*puntos, fill=color, width=2, smooth=False)
            canvas.create_text(margen-5, margen, text=f"{y_max:.3g}", anchor="e", fill=self.colores["texto_secundario"])
            canvas.create_text(margen-5, alto-margen, text=f"{y_min:.3g}", anchor="e", fill=self.colores["texto_secundario"])
            canvas.create_text(margen, alto-12, text=f"{t_min:.2f}", anchor="w", fill=self.colores["texto_secundario"])
            canvas.create_text(ancho-margen, alto-12, text=f"{t_max:.2f} t", anchor="e", fill=self.colores["texto_secundario"])

        canvas.bind("<Configure>", dibujar)
        self.after_idle(dibujar)

    def mostrar_analisis(self):
        if not self.modelo:
            self.mostrar_captura()
            return
        self.limpiar_pantalla()
        self.crear_encabezado(
            "Análisis de sensibilidad y trayectorias",
            "Compara cuánto cambian las estimaciones al hacer pasos de tiempo más pequeños.",
        )
        panel = tk.Frame(self, bg=self.colores["bg_tarjeta"], bd=1, relief="solid")
        panel.place(relx=0.04, rely=0.17, relwidth=0.92, relheight=0.68)
        tabla = tk.Frame(panel, bg=self.colores["bg_tarjeta"])
        tabla.pack(fill="x", padx=12, pady=(8, 2))
        filas = [
            (
                f"{dato['h']:.2f}", f"{dato['Y_final']:.6f}",
                f"{dato['theta_final']:.6f}", f"{dato['S']:.7f}",
                "—" if dato["diferencia_S"] is None else f"{dato['diferencia_S']:.3g}",
            )
            for dato in self.modelo["sensibilidad"]["pasos"]
        ]
        self.crear_tabla_resultados(
            tabla, filas, (
                "Tamaño del paso", "Respuesta final", "Temperatura final",
                "Índice general", "Cambio del índice",
            )
        )
        texto_tolerancias = (
            "El tamaño del paso es el intervalo entre cálculos: reducirlo permite "
            "observar si los resultados se mantienen parecidos. "
            "La diferencia del índice compara cada fila con el paso siguiente más pequeño."
            "\n\nTambién se probó la búsqueda de C con distintos márgenes de error: "
            + " | ".join(
            f"margen {dato['tolerancia']:.0e}: solución {dato['C']:.7f}, "
            f"{dato['iteraciones']} intentos, error restante "
            f"{abs(dato['residuo']):.2g}, índice general {dato['S']:.6f}"
            for dato in self.modelo["sensibilidad"]["tolerancias"]
            )
        )
        tk.Label(
            panel, text=texto_tolerancias, bg=self.colores["bg_tarjeta"],
            fg=self.colores["texto_secundario"], font=("Segoe UI", 9),
            wraplength=950, justify="left",
        ).pack(fill="x", padx=18, pady=(2, 4))
        graficas = tk.Frame(panel, bg=self.colores["bg_tarjeta"])
        graficas.pack(fill="both", expand=True, padx=12, pady=8)
        grafica_y = tk.Canvas(graficas, bg=self.colores["bg_tarjeta"], highlightthickness=0)
        grafica_theta = tk.Canvas(graficas, bg=self.colores["bg_tarjeta"], highlightthickness=0)
        grafica_y.pack(side="left", fill="both", expand=True, padx=4)
        grafica_theta.pack(side="left", fill="both", expand=True, padx=4)
        self.dibujar_grafica(
            grafica_y, self.modelo["documento_e"]["historial"],
            "Documento E: Y(t)", self.colores["primario"]
        )
        self.dibujar_grafica(
            grafica_theta, self.modelo["documento_f"]["historial"],
            "Documento F: temperatura estimada", "#f97316"
        )
        tk.Label(
            panel,
            text="En las gráficas, el eje horizontal representa el tiempo y el vertical muestra la respuesta o temperatura estimada.",
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
        self.configure(bg=self.colores["bg_principal"])
        self.crear_encabezado(
            "Operaciones e iteraciones detalladas",
            "Resultados intermedios, criterios de paro, supuestos y comparación de pasos."
        )
        contenedor = tk.Frame(self, bg=self.colores["bg_tarjeta"], bd=1, relief="solid")
        contenedor.place(relx=0.03, rely=0.16, relwidth=0.94, relheight=0.72)
        texto_area = tk.Text(
            contenedor, bg=self.colores["bg_tarjeta"], fg=self.colores["texto_oscuro"],
            font=("Consolas", 10), relief="flat", padx=20, pady=20, wrap="word",
            insertbackground="white",
        )
        scrollbar = tk.Scrollbar(contenedor, command=texto_area.yview)
        texto_area.configure(yscrollcommand=scrollbar.set)
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
            f"C converge cuando |x_nueva-x|<1e-10, máximo 100 iteraciones; "
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
        lineas.append("Ar se obtiene por Newton-Raphson desde w0=0.1; la raíz es la menor positiva.")
        lineas.append("Criterio para Ar: |w_nueva-w|<1e-12, máximo 99 iteraciones.")
        for paso in coefs_b["historial_ar"]:
            lineas.append(
                f"  iter {paso['iter']}: w={paso['x']:.10g}, "
                f"f(w)={paso['fx']:.4g}, f'(w)={paso['dfx']:.4g}, "
                f"w_nueva={paso['x_nueva']:.10g}"
            )
        lineas.append(
            "Cr usa mínimos cuadrados sobre (1,4.2),(2,6.1),(4,9.8),(6,13.9); "
            "Dr interpola (1,2.8),(2,4.6),(3,7.9),(4,12.5) en 2.5."
        )
        lineas.extend(("", "DOCUMENTO C — saturación"))
        lineas.extend(
            f"{clave}={val_c[clave]:.10g}" for clave in ("a_s", "b_s", "c_s", "d_s")
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
        lineas.append("Euler para p: y'=-0.6y+4t, y(0)=1, h=0.1:")
        lineas.extend(
            f"  t={tiempo:.1f}, y={valor:.10g}"
            for tiempo, valor in datos_d["historial_p"]
        )
        lineas.extend(("", "DOCUMENTO E — EDO de rendimiento"))
        lineas.append("Euler hacia adelante, intervalo [0,3], h=0.1; kappa se integra con Simpson 3/8.")
        lineas.append(
            f"r_sigma={val_e['r_sigma']}, sigma={val_e['sigma']}, "
            f"tau={val_e['tau']:.8g}, kappa={val_e['kappa']:.8g}, "
            f"omega={val_e['omega']:.8g}, h={val_e['h']}"
        )
        lineas.append(f"Y(0)={val_e['Y0']:.8g}; Y(3)={val_e['Y_final']:.8g}")
        lineas.extend(("", "DOCUMENTO F — EDO térmica"))
        lineas.append("Euler hacia adelante, intervalo [0,3], h=0.1; la trayectoria Y(t) de E se interpola linealmente entre pasos.")
        lineas.append(
            f"chi={val_f['chi']:.8g}, psi={val_f['psi']}, xi=a+d+g={val_f['xi']:.8g}"
        )
        lineas.append(f"theta(0)={val_f['theta0']} °C; theta(3)={val_f['theta_final']:.8g} °C")
        lineas.append(f"Justificación: {val_f['justificacion_theta0']}")
        lineas.append(
            f"Normalización λ={val_f['escala_fuente']:.8g}: "
            f"máx(Y/xi)={val_f['fuente_maxima']:.8g}; "
            f"capacidad máxima chi*psi²/4={val_f['capacidad_maxima']:.8g}; "
            f"la fuente se limita al {100*val_f['fraccion_capacidad']:.0f}% de esa capacidad."
        )
        lineas.extend(("", "DOCUMENTO G — índice global"))
        lineas.append(f"Raíces reales de saturación: {val_g['raices_saturacion']}")
        lineas.append(
            f"Se elige U*={val_g['U_star']:.10g}; integral Y(0,2)="
            f"{val_g['integral_Y']:.10g}; theta(1)={val_g['theta_1']:.10g}; S={val_g['S']:.10g}"
        )
        lineas.extend(("", "ANÁLISIS DE SENSIBILIDAD (Euler)"))
        for dato in self.modelo["sensibilidad"]["pasos"]:
            diferencia = "—" if dato["diferencia_S"] is None else f"{dato['diferencia_S']:.4g}"
            lineas.append(
                f"h={dato['h']:.3g}: Y(3)={dato['Y_final']:.8g}, "
                f"theta(3)={dato['theta_final']:.8g}, S={dato['S']:.8g}, "
                f"|ΔS|={diferencia}"
            )
        lineas.append("Sensibilidad a tolerancia de Newton-Raphson para C:")
        for dato in self.modelo["sensibilidad"]["tolerancias"]:
            lineas.append(
                f"  tol={dato['tolerancia']:.0e}: C={dato['C']:.10g}, "
                f"iteraciones={dato['iteraciones']}, residuo={dato['residuo']:.4g}, "
                f"S={dato['S']:.10g}, |ΔS|={dato['diferencia_S']:.4g}"
            )
        lineas.extend((
            "",
            "INTERPRETACIÓN",
            "S es un índice compuesto del expediente. Como no se proporcionan umbrales, no se etiqueta por sí solo como estable o inestable.",
            "",
            "SUPUESTOS DECLARADOS",
            "• Se usa r_sigma=1.2, la menor raíz positiva, para el decaimiento más conservador.",
            "• El expediente no define omega; se toma omega=kappa=L(3), por lo que kappa/omega=1.",
            "• Se propone theta(0)=25 °C como línea base ambiente/de reposo previa a la carga.",
            "• Se detectó que el término Y/xi original supera la capacidad logística chi*psi²/4 y hace divergir theta. Se añadió λ=max(1,max(Y/xi)/(0.8*chi*psi²/4)); el margen de 20% es una normalización de cálculo no dada en el PDF y debe describirse como tal.",
            "• Para U* se elige la menor raíz real positiva de la ecuación de saturación.",
        ))
        texto_area.insert("1.0", "\n".join(lineas))
        texto_area.configure(state="disabled")
        self.agregar_navegacion()


if __name__ == "__main__":
    Aplicacion().mainloop()