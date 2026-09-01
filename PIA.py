# Programa para capturar 4 matrículas de 7 dígitos,
# validar que no sean repetidas y que solo contengan números,
# y luego calcular N_i, P_i y S_i.

# Diccionario para guardar las matrículas como listas de dígitos
# Ejemplo: {"M1": [1, 9, 5, 0, 1, 9, 9], ...}
matriculas = {}

# Conjunto para evitar matrículas repetidas
matriculas_registradas = set()

# Bucle para capturar 4 matrículas
for i in range(1, 5):
    while True:
        valor = input(f"Ingrese la matrícula M{i} (7 dígitos): ")

        # Validación 1: solo números
        if not valor.isdigit() or len(valor) != 7:
            print("Error: la matrícula debe contener exactamente 7 dígitos numéricos.")
            continue

        # Validación 2: no repetir matrícula
        if valor in matriculas_registradas:
            print("Error: esa matrícula ya fue ingresada. Intenta con otra.")
            continue

        # Si pasa validaciones, se guarda
        matriculas_registradas.add(valor)
        matriculas[f"M{i}"] = [int(d) for d in valor]
        break

# Diccionario para guardar el valor completo de cada matrícula
# Ejemplo: {"M1": 1950199, "M2": 1234567, ...}
M = {
    clave: int("".join(str(d) for d in digitos))
    for clave, digitos in matriculas.items()
}

# Cálculo de N_i = M_i / 10^7
potencia_10_7 = 10 ** 7
N = {
    clave: valor / potencia_10_7
    for clave, valor in M.items()
}

# Cálculo de P_i = (1/7) * sum(digitos)
P = {
    clave: (1 / 7) * sum(digitos)
    for clave, digitos in matriculas.items()
}

# Cálculo de S_i = sum(k * dik) donde k es la posición del dígito
S = {
    clave: sum(pos * digito for pos, digito in enumerate(digitos, start=1))
    for clave, digitos in matriculas.items()
}

# Mostrar resultados
for clave in ["M1", "M2", "M3", "M4"]:
    print(f"{clave} = {matriculas[clave]}")
    print(f"Valor completo: {M[clave]}")
    print(f"N{clave[-1]} = {N[clave]}")
    print(f"P{clave[-1]} = {P[clave]}")
    print(f"S{clave[-1]} = {S[clave]}")
    print()