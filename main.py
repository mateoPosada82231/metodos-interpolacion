from sympy import symbols, expand, simplify, lambdify, Poly

def crear_funcion_lagrange(pares_xy):
    """
    Crea una función polinomial usando interpolación de Lagrange
    que pasa exactamente por todos los puntos dados

    Parámetros:
    pares_xy: lista de tuplas [(x1,y1), (x2,y2), ..., (xn,yn)]

    Retorna:
    - funcion: función de Python f(x) evaluable en cualquier punto
    - polinomio: expresión simbólica del polinomio
    - coeficientes: lista [a0, a1, a2, ...] para P(x) = a0 + a1*x + a2*x² + ...
    """

    x_vals = [p[0] for p in pares_xy]
    y_vals = [p[1] for p in pares_xy]
    x = symbols('x')
    polinomio = 0
    n = len(pares_xy)

    # Construir polinomio de Lagrange: P(x) = Σ yi * Li(x)
    for i in range(n):
        # Calcular Li(x) - polinomio base de Lagrange
        Li = 1
        for j in range(n):
            if i != j:
                Li *= (x - x_vals[j]) / (x_vals[i] - x_vals[j])
        polinomio += y_vals[i] * Li

    # Simplificar el polinomio
    polinomio_final = simplify(expand(polinomio))

    # Crear función evaluable
    funcion = lambdify(x, polinomio_final, 'numpy')

    # Extraer coeficientes
    try:
        poly_obj = Poly(polinomio_final, x)
        coeficientes = [float(c) for c in poly_obj.all_coeffs()[::-1]]
    except:
        coeficientes = [float(polinomio_final)] if polinomio_final.is_number else [0]

    return funcion, polinomio_final, coeficientes


def crear_funcion_newton(pares_xy):
    """
    Crea una función polinomial usando interpolación de Newton
    que pasa exactamente por todos los puntos dados

    Parámetros:
    pares_xy: lista de tuplas [(x1,y1), (x2,y2), ..., (xn,yn)]

    Retorna:
    - funcion: función de Python f(x) evaluable en cualquier punto
    - polinomio: expresión simbólica del polinomio
    - coeficientes: lista [a0, a1, a2, ...] para P(x) = a0 + a1*x + a2*x² + ...
    - tabla_diferencias: tabla de diferencias divididas
    """
    x_vals = [p[0] for p in pares_xy]
    y_vals = [p[1] for p in pares_xy]
    x = symbols('x')
    n = len(pares_xy)

    # Construir tabla de diferencias divididas
    tabla = [[0 for j in range(n)] for i in range(n)]

    # Primera columna - valores f(xi)
    for i in range(n):
        tabla[i][0] = y_vals[i]

    # Calcular diferencias divididas
    for j in range(1, n):
        for i in range(n - j):
            tabla[i][j] = (tabla[i + 1][j - 1] - tabla[i][j - 1]) / (x_vals[i + j] - x_vals[i])

    # Construir polinomio de Newton
    polinomio = tabla[0][0]  # f[x0]
    producto = 1

    for i in range(1, n):
        producto *= (x - x_vals[i - 1])
        polinomio += tabla[0][i] * producto

    # Simplificar
    polinomio_final = simplify(expand(polinomio))

    # Crear función evaluable
    funcion = lambdify(x, polinomio_final, 'numpy')

    # Extraer coeficientes
    try:
        poly_obj = Poly(polinomio_final, x)
        coeficientes = [float(c) for c in poly_obj.all_coeffs()[::-1]]
    except:
        coeficientes = [float(polinomio_final)] if polinomio_final.is_number else [0]

    return funcion, polinomio_final, coeficientes, tabla


def leer_puntos_csv(ruta):
    """
    Lee un archivo CSV con dos columnas separadas por espacio, usando coma como separador decimal.
    Devuelve una lista de tuplas (x, y) ordenadas por x.
    """
    puntos = []
    with open(ruta, 'r', encoding='utf-8') as f:
        for linea in f:
            if linea.strip():
                partes = linea.strip().split()
                if len(partes) == 2:
                    x = float(partes[0].replace(',', '.'))
                    y = float(partes[1].replace(',', '.'))
                    puntos.append((x, y))
    puntos.sort(key=lambda t: t[0])
    return puntos


# Definir los puntos por los que debe pasar la función
puntos = leer_puntos_csv('Datos/DatosDePrueba.csv')

# MÉTODO DE LAGRANGE - Crear función interpolante
f_lagrange, polinomio_lag, coef_lag = crear_funcion_lagrange(puntos)

print("🔹 LAGRANGE")
print(f"Polinomio: P(x) = {polinomio_lag.evalf(5)}")
print(f"Coeficientes: {[round(c, 5) for c in coef_lag]}")

# MÉTODO DE NEWTON - Crear función interpolante
f_newton, polinomio_new, coef_new, tabla = crear_funcion_newton(puntos)

print("🔹 NEWTON")
print(f"Polinomio: P(x) = {polinomio_new.evalf(5)}")
print(f"Coeficientes: {[round(c, 5) for c in coef_new]}")

# USAR LAS FUNCIONES CREADAS
print(f"\nEvaluaciones:")
print(f"f_lagrange(1) = {round(f_lagrange(1), 5)}")
print(f"f_newton(1) = {round(f_newton(1), 5)}")

# VERIFICAR que pasan por los puntos originales
print(f"\nVerificación:")
for x, y in puntos:
    resultado = round(f_lagrange(x), 5)
    print(f"Punto ({x}, {y}): f({x}) = {resultado} ✓")
