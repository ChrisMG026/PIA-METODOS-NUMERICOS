# Variables para guardar el valor total de cada matrícula
M1 = 0
M2 = 0
M3 = 0
M4 = 0

# Variables para guardar cada dígito de cada matrícula
di1_1 = 0
di2_1 = 0
di3_1 = 0
di4_1 = 0
di5_1 = 0
di6_1 = 0
di7_1 = 0

di1_2 = 0
di2_2 = 0
di3_2 = 0
di4_2 = 0
di5_2 = 0
di6_2 = 0
di7_2 = 0

di1_3 = 0
di2_3 = 0
di3_3 = 0
di4_3 = 0
di5_3 = 0
di6_3 = 0
di7_3 = 0

di1_4 = 0
di2_4 = 0
di3_4 = 0
di4_4 = 0
di5_4 = 0
di6_4 = 0
di7_4 = 0

# Potencia que se usará para dividir
potencia_10_7 = 10 ** 7

# Conjunto para evitar matrículas repetidas
matriculas_registradas = set()

for i in range(1, 5):
    while True:
        valor = input(f"Ingrese la matrícula M{i} (7 dígitos): ")
        
        # Validar que sea numérico y tenga 7 dígitos
        if not valor.isdigit() or len(valor) != 7:
            print("Error: la matrícula debe contener exactamente 7 dígitos numéricos.")
            continue
        
        # Validar que no se repita la matrícula
        if valor in matriculas_registradas:
            print("Error: esa matrícula ya fue ingresada. Intenta con otra.")
            continue
        
        # Si pasa las validaciones, guardar matrícula
        matriculas_registradas.add(valor)
        
        # Valor total de la matrícula
        matricula = int(valor)
        
        # Guardar el valor total de la matrícula
        globals()[f"M{i}"] = matricula
        
        # Guardar cada dígito como di1, di2, ..., di7
        digitos = [int(d) for d in valor]
        
        for j, digito in enumerate(digitos, start=1):
            globals()[f"di{j}_{i}"] = digito
        
        break

# Mostrar las matrículas completas
print(f"M1 = {M1}")
print(f"M2 = {M2}")
print(f"M3 = {M3}")
print(f"M4 = {M4}")

# Mostrar dígitos por matrícula
print(f"M1 = ({di1_1}, {di2_1}, {di3_1}, {di4_1}, {di5_1}, {di6_1}, {di7_1})")
print(f"M2 = ({di1_2}, {di2_2}, {di3_2}, {di4_2}, {di5_2}, {di6_2}, {di7_2})")
print(f"M3 = ({di1_3}, {di2_3}, {di3_3}, {di4_3}, {di5_3}, {di6_3}, {di7_3})")
print(f"M4 = ({di1_4}, {di2_4}, {di3_4}, {di4_4}, {di5_4}, {di6_4}, {di7_4})")

# Cálculo de N1, N2, N3, N4
N1 = M1 / potencia_10_7
N2 = M2 / potencia_10_7
N3 = M3 / potencia_10_7
N4 = M4 / potencia_10_7

# Mostrar resultados
print(f"N1 = {N1}")
print(f"N2 = {N2}")
print(f"N3 = {N3}")
print(f"N4 = {N4}")

# Cálculo de Si = sum(k * dik) para cada matrícula
S1 = (1 * di1_1) + (2 * di2_1) + (3 * di3_1) + (4 * di4_1) + (5 * di5_1) + (6 * di6_1) + (7 * di7_1)
S2 = (1 * di1_2) + (2 * di2_2) + (3 * di3_2) + (4 * di4_2) + (5 * di5_2) + (6 * di6_2) + (7 * di7_2)
S3 = (1 * di1_3) + (2 * di2_3) + (3 * di3_3) + (4 * di4_3) + (5 * di5_3) + (6 * di6_3) + (7 * di7_3)
S4 = (1 * di1_4) + (2 * di2_4) + (3 * di3_4) + (4 * di4_4) + (5 * di5_4) + (6 * di6_4) + (7 * di7_4)

print(f"S1 = {S1}")
print(f"S2 = {S2}")
print(f"S3 = {S3}")
print(f"S4 = {S4}")