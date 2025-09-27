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
    polinomio_sin_simplificar = polinomio
    polinomio_final = simplify(expand(polinomio))

    # Crear función evaluable
    funcion = lambdify(x, polinomio_final, 'numpy')

    # Extraer coeficientes
    try:
        poly_obj = Poly(polinomio_final, x)
        coeficientes = [float(c) for c in poly_obj.all_coeffs()[::-1]]
    except:
        coeficientes = [float(polinomio_final)] if polinomio_final.is_number else [0]

    return funcion, polinomio_final, coeficientes, polinomio_sin_simplificar


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
    polinomio_sin_simplificar = polinomio
    polinomio_final = simplify(expand(polinomio))

    # Crear función evaluable
    funcion = lambdify(x, polinomio_final, 'numpy')

    # Extraer coeficientes
    try:
        poly_obj = Poly(polinomio_final, x)
        coeficientes = [float(c) for c in poly_obj.all_coeffs()[::-1]]
    except:
        coeficientes = [float(polinomio_final)] if polinomio_final.is_number else [0]

    return funcion, polinomio_final, coeficientes, tabla, polinomio_sin_simplificar

def escoger_ruta():
    """
    Solicita al usuario la ruta del archivo CSV y la devuelve.
    """
    des = ""
    while des not in ["1Y", "1L", "2Y", "2L", "3Y", "3L"]:
        des = input("Elige el conjunto de datos (1Y, 1L, 2Y, 2L, 3Y, 3L): ").strip().upper()
        if des not in ["1Y", "1L", "2Y", "2L", "3Y", "3L"]:
            print("Opción no válida. Intenta de nuevo.")
    ruta = f"Datos/{des}.csv"
    return ruta

# python
def leer_puntos_csv(ruta):
    """
    Lee un archivo CSV con dos columnas separadas por coma.
    Devuelve una lista de tuplas (x, y) ordenadas por x, con decimales limitados a 5.
    """
    puntos = []
    with open(ruta, 'r', encoding='utf-8') as f:
        for linea in f:
            if linea.strip():
                partes = linea.strip().split(',')
                if len(partes) == 2:
                    x_str = partes[0].replace(',', '.').rstrip('.')
                    y_str = partes[1].replace(',', '.').rstrip('.')
                    x = round(float(x_str), 5)
                    y = round(float(y_str), 5)
                    puntos.append((x, y))
    puntos.sort(key=lambda t: t[0])
    return puntos


def procesar_y_reportar(ruta_archivo, stream):
    """
    Realiza el análisis de interpolación para un archivo de datos y escribe
    los resultados en un stream (puede ser la consola o un archivo).
    """
    stream.write(f"\n{'='*80}\n")
    stream.write(f"REPORTE PARA EL ARCHIVO: {ruta_archivo}\n")
    stream.write(f"{'='*80}\n")

    todos_los_puntos = leer_puntos_csv(ruta_archivo)

    # Seleccionar puntos específicos para la interpolación (índices 1-based)
    indices_a_tomar = [7, 10, 13, 16, 19, 21, 24]
    puntos_interpolacion = [todos_los_puntos[i - 1] for i in indices_a_tomar if i <= len(todos_los_puntos)]

    # Puntos restantes para la tabla de errores
    puntos_verificacion = [p for i, p in enumerate(todos_los_puntos) if (i + 1) not in indices_a_tomar]

    stream.write(f"\nSe leyeron {len(todos_los_puntos)} puntos de {ruta_archivo}.\n")
    stream.write(f"Se usarán {len(puntos_interpolacion)} puntos para la interpolación.\n")

    # Mostrar puntos tomados y el intervalo
    stream.write("\nPuntos utilizados para la interpolación:\n")
    stream.write(str(puntos_interpolacion) + "\n")
    intervalo_x = (puntos_interpolacion[0][0], puntos_interpolacion[-1][0])
    stream.write(f"Intervalo en X de los puntos de interpolación: {intervalo_x}\n")

    # MÉTODO DE LAGRANGE
    f_lagrange, polinomio_lag, coef_lag, pol_lag_sin_simplificar = crear_funcion_lagrange(puntos_interpolacion)
    stream.write("\n🔹 LAGRANGE\n")
    stream.write(f"Polinomio sin simplificar: P(x) = {pol_lag_sin_simplificar.evalf(5)}\n")
    stream.write(f"Polinomio simplificado: P(x) = {polinomio_lag.evalf(5)}\n")
    stream.write(f"Coeficientes: {[round(c, 5) for c in coef_lag]}\n")

    # MÉTODO DE NEWTON
    f_newton, polinomio_new, coef_new, tabla, pol_new_sin_simplificar = crear_funcion_newton(puntos_interpolacion)
    stream.write("\n🔹 NEWTON\n")
    stream.write(f"Polinomio sin simplificar: P(x) = {pol_new_sin_simplificar.evalf(5)}\n")
    stream.write(f"Polinomio simplificado: P(x) = {polinomio_new.evalf(5)}\n")
    stream.write(f"Coeficientes: {[round(c, 5) for c in coef_new]}\n")

    # TABLA DE ERRORES
    stream.write("\n🔹 TABLA DE ERRORES (puntos no utilizados en la interpolación)\n")
    stream.write("-" * 55 + "\n")
    stream.write(f"{'Punto (x, y)':<25} | {'Valor Polinomio P(x)':<20} | {'Error |y - P(x)|':<15}\n")
    stream.write("-" * 55 + "\n")
    for x, y in puntos_verificacion:
        valor_evaluado = f_lagrange(x)
        error = abs(y - valor_evaluado)
        stream.write(f"({x:<10}, {y:<10}) | {round(valor_evaluado, 5):<20} | {round(error, 5):<15}\n")
    stream.write("-" * 55 + "\n")

    # VERIFICACIÓN
    stream.write(f"\nVerificación (puntos de interpolación):\n")
    for x, y in puntos_interpolacion:
        resultado = round(f_lagrange(x), 5)
        stream.write(f"Punto ({x}, {y}): f({x}) = {resultado} (Error: {round(abs(y - resultado), 5)}) ✓\n")
    stream.write("\n")


if __name__ == "__main__":
    # Lista de todos los archivos de datos a procesar
    archivos_datos = [
        "Datos/1Y.csv", "Datos/1L.csv",
        "Datos/2Y.csv", "Datos/2L.csv",
        "Datos/3Y.csv", "Datos/3L.csv"
    ]

    # Abrir el archivo de reporte en modo escritura
    with open("reporte.txt", "w", encoding="utf-8") as f_reporte:
        # Procesar cada archivo y escribir en el reporte
        for archivo in archivos_datos:
            try:
                procesar_y_reportar(archivo, f_reporte)
            except Exception as e:
                f_reporte.write(f"\n{'!'*80}\n")
                f_reporte.write(f"ERROR al procesar el archivo {archivo}: {e}\n")
                f_reporte.write(f"{'!'*80}\n\n")

    print("El reporte 'reporte.txt' ha sido generado exitosamente.")
