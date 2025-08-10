import psutil
import os
import inspect
import time
import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from contextlib import redirect_stdout

def medir_funcao(func, *args, **kwargs):
    p = psutil.Process(os.getpid())

    cpu_times_start = p.cpu_times()
    mem_start = p.memory_info()

    inicio = time.perf_counter()
    resultado = func(*args, **kwargs)
    fim = time.perf_counter()

    cpu_times_end = p.cpu_times()
    mem_end = p.memory_info()

    mem_diff_bytes = mem_end.rss - mem_start.rss
    mem_diff_kb = mem_diff_bytes / 1024

    print(f"Tempo real: {fim - inicio:.6f} s")
    print(f"Tempo de CPU (user): {cpu_times_end.user - cpu_times_start.user:.20f} s")
    print(f"Memória usada: {mem_diff_bytes} bytes ({mem_diff_kb:.2f} KB)")

    return resultado


def resolver_edo(expr_edo, condicoes_iniciais, x_val, n, h):
    x = sp.Symbol('x')
    y = sp.Function('y')

    solucao = sp.dsolve(expr_edo, ics=condicoes_iniciais)

    x_vals = [x_val + i * h for i in range(n + 1)]
    y_vals = [solucao.rhs.subs(x, val).evalf() for val in x_vals]

    return x_vals, y_vals

# Euler
# Vantagens: Baixo Custo
# Desvantagens: Precisa de passos pequenos para ter um erro pequeno. Pode se tornar instável se o passo for muito grande
# Caso de uso: fins didáticos
def euler(f, y0, t0, h, n):
    y, t = [y0], [t0]
    for _ in range(n):
        y0 += h * f(t0, y0)
        t0 += h
        y.append(y0)
        t.append(t0)
    return t, y

# Euler Implícito
# Vantagens: É mais estável em passos maiores. Por conta da estabilidade, pode ser usado em longos períodos
# Desvantagens: Custo computacional alto
# Caso de uso: EDO rígido
def euler_implicito(f, y0, t0, h, n):
    from scipy.optimize import fsolve
    y, t = [y0], [t0]
    for _ in range(n):
        t1 = t0 + h
        g = lambda y1: y1 - y0 - h * float(f(float(t1), float(y1)))
        y1 = fsolve(g, y0)[0]
        y0, t0 = y1, t1
        y.append(y0)
        t.append(t0)
    return t, y


# RK2 - Ponto Médio
# Vantagens: Precisão melhor do que Euler. Não exige a resolução de equações implícitas
# Desvantagens: Custo alto. Não adequado para EDOs Rígidas
# Caso de uso: Problemas não rígidos com precisão moderado
def rk2_ponto_medio(f, y0, t0, h, n):
    y, t = [y0], [t0]
    for _ in range(n):
        k1 = f(t0, y0)
        k2 = f(t0 + h / 2, y0 + h / 2 * k1)
        y0 += h * k2
        t0 += h
        y.append(y0)
        t.append(t0)
    return t, y

# RK2 - Heun
# Vantagens: Melhor precisão que Euler, Explícito
# Desvantagens: Custo alto. Não adequado para EDO rígida
# Caso de uso: Problemas não rígidos com precisão moderada.
def rk2_heun(f, y0, t0, h, n):
    y, t = [y0], [t0]
    for _ in range(n):
        k1 = f(t0, y0)
        k2 = f(t0 + h, y0 + h * k1)
        y0 += h / 2 * (k1 + k2)
        t0 += h
        y.append(y0)
        t.append(t0)
    return t, y

# Trapezio Implícito
# Vantagens: É estável
# Desvantagens: Custo computacional por passo
# Caso de uso: EDO rígidas com boa precisão
def trapezio_implicito(f, y0, t0, h, n):
    from scipy.optimize import fsolve
    y, t = [y0], [t0]
    for _ in range(n):
        t1 = t0 + h
        f0 = f(float(t0), float(y0))  # garantir float
        g = lambda y1: y1 - y0 - h/2 * (f0 + f(float(t1), float(y1)))  # garantir float em f
        y1 = fsolve(g, y0)[0]
        y0, t0 = y1, t1
        y.append(y0)
        t.append(t0)
    return t, y

# BDF2 Implícito
# Vantagens: Excelente estabilidade para EDO rígidas.
# Desvantagens: Custo Alto, necessita de pontos inicias
# Caso de uso: EDOs Rígidas de longos períodos
def bdf2(f, y0, y1, t0, h, n):
    from scipy.optimize import fsolve
    y = [y0, y1]
    t = [t0, t0 + h]
    for k in range(1, n):
        t_next = t[k] + h
        g = lambda y_next: y_next - (4/3) * y[k] + (1/3) * y[k-1] - (2/3)*h * float(f(t_next, y_next))
        y_next = fsolve(g, y[k])[0]
        y.append(y_next)
        t.append(t_next)
    return t, y


# Adams-Bashfort - segunda

# Vantagens: Rápido, baixo custo
# Desvantagens: Instável com passo grande, necessita de um passo anterior
# Caso de uso: Simulações de curto prazo ou onde o custo computacional é mais importante que precisão ou robustez.
def adams_bashforth_2(f, y0, t0, h, n):
    y = [y0]
    t = [t0]
    y.append(y0 + h * f(t0, y0))
    t.append(t0 + h)

    for k in range(1, n):
        fk = f(t[k], y[k])
        fk_1 = f(t[k-1], y[k-1])
        y_next = y[k] + h * (3/2 * fk - 1/2 * fk_1)
        y.append(y_next)
        t.append(t[k] + h)
    return t, y

# Adams-Moulton - segunda
# Vantagens: Baixo Custo, mais preciso que o segunda ordem
# Desvantagens: Muito sensível a instabilidade, Necessita de 2 valores prévios
# Caso de uso: Quando se deseja mais precisão que o de segunda ordem
def adams_moulton_2(f, y0, t0, h, n):
    from scipy.optimize import fsolve
    y = [y0]
    t = [t0]
    y.append(y0 + h * f(t0, y0))
    t.append(t0 + h)

    for k in range(1, n):
        t_next = t[k] + h
        fk = f(t[k], y[k])
        fk_1 = f(t[k-1], y[k-1])
        g = lambda y_next: y[k] + h * (5/12 * f(t_next, y_next) + 2/3 * fk - 1/12 * fk_1) - y_next
        y_next = fsolve(g, y[k])[0]
        y.append(y_next)
        t.append(t_next)
    return t, y

# Adams-Bashforth - terceira
# Vantagens: Alta precisão
# Desvantagens: Instavel em passos largos, Necessita de 4 valores prévios
# Caso de uso: Simulações onde o problema é suave e você quer velocidade com boa precisão.
def adams_bashforth_3(f, y0, t0, h, n):
    y = [y0]
    t = [t0]
    y.append(y0 + h * f(t0, y0))
    t.append(t0 + h)
    y.append(y[1] + h * f(t[1], y[1]))
    t.append(t[1] + h)

    for k in range(2, n):
        fk = f(t[k], y[k])
        fk_1 = f(t[k-1], y[k-1])
        fk_2 = f(t[k-2], y[k-2])
        y_next = y[k] + h * (23/12 * fk - 4/3 * fk_1 + 5/12 * fk_2)
        y.append(y_next)
        t.append(t[k] + h)
    return t, y

# RK3
# Vantagens: Simples e direto
# Desvantagens: Custo moderado. Instável em problemas stiff.
# Caso de uso: Ideal para problemas não rígidos, onde simplicidade e desempenho são mais importantes que máxima precisão. Ótimo para Inicialização.
def rk3(f, y0, t0, h, n):
    y = [y0]
    t = [t0]
    for _ in range(n):
        F1 = f(t0, y0)
        F2 = f(t0 + h/2, y0 + h/2 * F1)
        F3 = f(t0 + h, y0 - h * F1 + 2 * h * F2)
        y0 += h / 6 * (F1 + 4 * F2 + F3)
        t0 += h
        y.append(y0)
        t.append(t0)
    return t, y

# Adams-Moulton - terceira
# Vantagens: alta precisão, funciona com passos maiores, iplícito e estável
# Desvantagens: Alto custo computacional. Pode sr sensível à escolha iniceial do chute
# Caso de uso: Ideal para problemas rígidos (reação química, circuitos elétricos, etc.). Quando se busca alta estabilidade e se pode tolerar maior custo por passo.
def adams_moulton_3(f, y0, t0, h, n):
    from scipy.optimize import fsolve
    y = [y0]
    t = [t0]
    y.append(y0 + h * f(t0, y0))
    t.append(t0 + h)
    y.append(y[1] + h * f(t[1], y[1]))
    t.append(t[1] + h)

    for k in range(2, n):
        t_next = t[k] + h
        fk = f(t[k], y[k])
        fk_1 = f(t[k-1], y[k-1])
        fk_2 = f(t[k-2], y[k-2])
        g = lambda y_next: y[k] + h * (3/8 * f(t_next, y_next) + 19/24 * fk - 5/24 * fk_1 + 1/24 * fk_2) - y_next
        y_next = fsolve(g, y[k])[0]
        y.append(y_next)
        t.append(t_next)
    return t, y

# Adams-Bashforth - quarta
# Vantagens: Baixo custo por passo e rápido
# Desvantagens: Instável para problemas rígidos ou com passos grandes. Erros numéricos podem se propagar fortemente.
# Caso de uso: Simulações de sistemas suaves, sem rigidez. alta velocidade e precisão moderada.
def adams_bashforth_4(f, y0, t0, h, n):
    y = [y0]
    t = [t0]
    for _ in range(3):
        y.append(y[-1] + h * f(t[-1], y[-1]))
        t.append(t[-1] + h)

    for k in range(3, n):
        fk = f(t[k], y[k])
        fk_1 = f(t[k-1], y[k-1])
        fk_2 = f(t[k-2], y[k-2])
        fk_3 = f(t[k-3], y[k-3])
        y_next = y[k] + h * (55/24 * fk - 59/24 * fk_1 + 37/24 * fk_2 - 3/8 * fk_3)
        y.append(y_next)
        t.append(t[k] + h)
    return t, y

def rk4(f, y0, t0, h, n):
    t = [t0]
    y = [y0]
    for _ in range(n):
        k1 = f(t[-1], y[-1])
        k2 = f(t[-1] + h/2, y[-1] + h*k1/2)
        k3 = f(t[-1] + h/2, y[-1] + h*k2/2)
        k4 = f(t[-1] + h, y[-1] + h*k3)
        y_next = y[-1] + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
        t_next = t[-1] + h
        t.append(t_next)
        y.append(y_next)
    return t, y


# Método preditor-corretor
# Vantagens: Alta Precisão, Estável e passo adaptável
# Desvantagens: Custo computacional, Acumulo de erro
# Casos de uso: Quando você quer mais precisão e estabilidade que AB4, mas não quer pagar o preço de resolver sistemas como no AM3. Situação com controle de custo computacional
def adams_bashforth_moulton_pc4(f, y0, t0, tf, h, num_corrector_iterations=1):
    num_initial_steps = 3
    t_init, y_init = rk4(f, y0, t0, h, num_initial_steps)

    t_values = list(t_init)
    y_values = list(y_init)

    f_values_hist = [f(t_values[i], y_values[i]) for i in range(num_initial_steps + 1)]

    current_idx = num_initial_steps

    while t_values[current_idx] < tf:
        t_n = t_values[current_idx]
        y_n = y_values[current_idx]
        t_next = t_n + h

        if t_next > tf:
            h_actual = tf - t_n
            if h_actual <= 0:
                break

            pass

        f_n_minus_3 = f_values_hist[0]
        f_n_minus_2 = f_values_hist[1]
        f_n_minus_1 = f_values_hist[2]
        f_n = f_values_hist[3]

        y_predicted = y_n + (h / 24.0) * (55 * f_n - 59 * f_n_minus_1 + 37 * f_n_minus_2 - 9 * f_n_minus_3)

        y_corrected = y_predicted

        for _ in range(num_corrector_iterations):
            f_predicted_or_corrected = f(t_next, y_corrected)
            y_corrected = y_n + (h / 24.0) * (9 * f_predicted_or_corrected + 19 * f_n - 5 * f_n_minus_1 + f_n_minus_2)

        y_final_step = y_corrected

        t_values.append(t_next)
        y_values.append(y_final_step)

        f_values_hist.pop(0)
        f_values_hist.append(f(t_next, y_final_step))

        current_idx += 1

    return t_values, y_values

# Dormand Prince
# Vantagens: Controle do erro de forma adaptativa, Alta precisão, Ele não precisa de métodos de inicialização
# Desvantagens: Custo computacional, Ineficiente em EDOs rígidas
# Casos de uso: Problemas não-rigidos com comportamento variável, Aplicações que precisa controlar o erro
def dormand_prince_54(f, y0, t0, tf, h_initial, atol=1e-7, rtol=1e-5):
    c2, c3, c4, c5, c6, c7 = 1/5, 3/10, 4/5, 8/9, 1, 1

    a21 = 1/5
    a31, a32 = 3/40, 9/40
    a41, a42, a43 = 44/45, -56/15, 32/9
    a51, a52, a53, a54 = 19372/6561, -25360/2187, 64448/6561, -212/729
    a61, a62, a63, a64, a65 = 9017/3168, -355/33, 46732/5247, 49/176, -5103/18656
    a71, a72, a73, a74, a75, a76 = 35/384, 0, 500/1113, 125/192, -2187/6784, 11/84

    b5_1, b5_2, b5_3, b5_4, b5_5, b5_6, b5_7 = 35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0

    b4_hat_1, b4_hat_2, b4_hat_3, b4_hat_4, b4_hat_5, b4_hat_6, b4_hat_7 = 5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40

    t_values = []
    y_values = []

    t = t0
    if isinstance(y0, (int, float)):
        y = np.array([y0])
    else:
        y = np.array(y0)

    h = h_initial

    k1 = f(t, y)

    t_values.append(t)
    y_values.append(y.copy())

    while t < tf:
        if t + h > tf:
            h = tf - t
            if h < 1e-12:
                break

        k2 = f(t + c2 * h, y + h * a21 * k1)
        k3 = f(t + c3 * h, y + h * (a31 * k1 + a32 * k2))
        k4 = f(t + c4 * h, y + h * (a41 * k1 + a42 * k2 + a43 * k3))
        k5 = f(t + c5 * h, y + h * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4))
        k6 = f(t + c6 * h, y + h * (a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5))
        k7 = f(t + c7 * h, y + h * (a71 * k1 + a72 * k2 + a73 * k3 + a74 * k4 + a75 * k5 + a76 * k6))

        y_next_5 = y + h * (b5_1*k1 + b5_2*k2 + b5_3*k3 + b5_4*k4 + b5_5*k5 + b5_6*k6 + b5_7*k7)
        y_next_4 = y + h * (b4_hat_1*k1 + b4_hat_2*k2 + b4_hat_3*k3 + b4_hat_4*k4 + b4_hat_5*k5 + b4_hat_6*k6 + b4_hat_7*k7)

        abs_error_estimate = np.max(np.abs(y_next_5 - y_next_4))

        if y.size > 1:
            tol = atol + np.linalg.norm(y_next_5) * rtol
        else:
            tol = atol + np.abs(y_next_5[0]) * rtol

        safety_factor = 0.9
        power = 1.0 / 5.0

        if abs_error_estimate <= 1e-18:
            s = 2.0
        else:
            s = safety_factor * (tol / abs_error_estimate)**power

        s = max(0.2, min(s, 10.0))

        if abs_error_estimate <= tol:
            t += h
            y = y_next_5
            k1 = k7

            t_values.append(t)
            y_values.append(y.copy())

            h = h * s
        else:
            h = h * s
            k1 = f(t, y)

    return t_values, y_values

# Runge-Kutta-Fehlberg
# Vantagens: Controle de erro em cada passo e ajusto-o conforme necessário, Alta precisão, Ele não precisa de métodos de inicialização
# Desvantagens: Custo computacional, Ineficiente em EDOs rígidas
# Casos de uso: Problemas não-rigidos com comportamento variável, Aplicações que precisa controlar o erro
def rkf45(f, y0, t0, tf, h_initial, atol=1e-7, rtol=1e-5):
    c = np.array([0, 1/4, 3/8, 12/13, 1, 1/2])

    a21 = 1/4
    a31, a32 = 3/32, 9/32
    a41, a42, a43 = 1932/2197, -7200/2197, 7296/2197
    a51, a52, a53, a54 = 439/216, -8, 3680/513, -845/4104
    a61, a62, a63, a64, a65 = -8/27, 2, -3544/2565, 1859/4104, -11/40

    b5 = np.array([16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55])
    b4_hat = np.array([25/216, 0, 1408/2565, 2197/4104, -1/5, 0])

    t_values = []
    y_values = []

    t = t0
    if isinstance(y0, (int, float)):
        y = np.array([y0])
    else:
        y = np.array(y0)

    h = h_initial

    t_values.append(t)
    y_values.append(y.copy())

    while t < tf:
        if t + h > tf:
            h = tf - t
            if h < 1e-12:
                break

        k1 = f(t, y)
        k2 = f(t + a21 * h, y + a21 * h * k1)
        k3 = f(t + a31 * h + a32 * h, y + a31 * h * k1 + a32 * h * k2)
        k4 = f(t + a41 * h + a42 * h + a43 * h, y + a41 * h * k1 + a42 * h * k2 + a43 * h * k3)
        k5 = f(t + a51 * h + a52 * h + a53 * h + a54 * h, y + a51 * h * k1 + a52 * h * k2 + a53 * h * k3 + a54 * h * k4)
        k6 = f(t + a61 * h + a62 * h + a63 * h + a64 * h + a65 * h, y + a61 * h * k1 + a62 * h * k2 + a63 * h * k3 + a64 * h * k4 + a65 * h * k5)

        y_next_5 = y + h * (b5[0]*k1 + b5[1]*k2 + b5[2]*k3 + b5[3]*k4 + b5[4]*k5 + b5[5]*k6)
        y_next_4 = y + h * (b4_hat[0]*k1 + b4_hat[1]*k2 + b4_hat[2]*k3 + b4_hat[3]*k4 + b4_hat[4]*k5 + b4_hat[5]*k6)

        error_estimate = np.linalg.norm(y_next_5 - y_next_4)

        if y.size > 1:
            tol = atol + np.linalg.norm(y_next_5) * rtol
        else:
            tol = atol + np.abs(y_next_5[0]) * rtol

        safety_factor = 0.9
        power = 1.0 / 5.0

        if error_estimate < 1e-18:
            s = 2.0
        else:
            s = safety_factor * (tol / error_estimate)**power

        s = max(0.2, min(s, 5.0))

        if error_estimate <= tol:
            t += h
            y = y_next_5

            t_values.append(t)
            y_values.append(y.copy())

            h = h * s
        else:
            h = h * s

    return t_values, y_values

x = sp.Symbol('x')
y = sp.Function('y')

def edo_para_funcao(expr_edo, y_sym, x_sym):
    deriv = sp.Derivative(y_sym, x_sym)
    f_expr = sp.solve(expr_edo, deriv)[0]

    def f_func(t, y):
        t_val = float(t)
        y_val = float(y)
        return float(f_expr.subs({x_sym: t_val, y_sym: y_val}))
    
    return f_func

def limpar_nome(nome):
            return nome.replace(" ", "_").replace("(", "").replace(")", "").replace("/", "_")

edos = [
    {
        "titulo": "EDO Crescimento Populacional com capacidade de carga",
        "expr_edo": sp.Eq(sp.Derivative(y(x), x) - 0.5 * y(x) * (1 - (y(x) / 1000)), 0),
        "condicoes_iniciais": {y(0): 100},
        "x_val": 12,
        "passo": 1,
    },
    {
        "titulo": "EDO de resfriamento de um objeto de Newton",
        "expr_edo": sp.Eq(sp.Derivative(y(x), x) + 0.1 * (y(x) - 20), 0),
        "condicoes_iniciais": {y(0): 150},
        "x_val": 10,
        "passo": 2,
    },
    {
        "titulo": "EDO Reação Química de Primeira Ordem",
        "expr_edo": sp.Eq(sp.Derivative(y(x), x) + 0.02 * y(x), 0),
        "condicoes_iniciais": {y(0): 1},
        "x_val": 20,
        "passo": 2,
    },
    {
        "titulo": "EDO Decaimento Radioativo",
        "expr_edo": sp.Eq(sp.Derivative(y(x), x) + 0.04 * y(x), 0),
        "condicoes_iniciais": {y(0): 500},
        "x_val": 5,
        "passo": 0.5,
    },
    {
        "titulo": "EDO Aquisição de Habilidades",
        "expr_edo": sp.Eq(sp.Derivative(y(x), x) - 0.2 * (10 - y(x)), 0),
        "condicoes_iniciais": {y(0): 1},
        "x_val": 8,
        "passo": 0.8,
    },
    {
        "titulo": "EDO de esfriamento de xícara de café em um ambiente variável",
        "expr_edo": sp.Eq(sp.Derivative(y(x), x) + 0.07 * (y(x) - (20 + 5 * sp.sin((sp.pi * y(x))/y(x)))), 0),
        "condicoes_iniciais": {y(0): 90},
        "x_val": 30,
        "passo": 2,
    },
    {
        "titulo": "EDO de corrente em Circuito RL Simples",
        "expr_edo": sp.Eq(sp.Derivative(y(x), x) - 10 + 20 * y(x), 0),
        "condicoes_iniciais": {y(0): 0},
        "x_val": 3,
        "passo": 0.3,
    },
    {
        "titulo": "EDO de Propagação de Doenças (SIR simplificado)",
        "expr_edo": sp.Eq(sp.Derivative(y(x), x) + 0.005 * y(x), 0),
        "condicoes_iniciais": {y(0): 500},
        "x_val": 10,
        "passo": 1,
    },
    {
        "titulo": "EDO de Crescimento de População de Peixes em Lago",
        "expr_edo": sp.Eq(sp.Derivative(y(x), x) - 0.3 * y(x) * sp.log(5000/ y(x)), 0),
        "condicoes_iniciais": {y(0): 100},
        "x_val": 5,
        "passo": 0.5,
    },
    {
        "titulo": "EDO de descarga de um capacitor em Circuito RC",
        "expr_edo": sp.Eq(sp.Derivative(y(x), x) - 10 * (12 * sp.cos(2 * sp.pi * x) - y(x)), 0),
        "condicoes_iniciais": {y(0): 0},
        "x_val": 5,
        "passo": 0.25,
    },
]

metodos = [
    ("Euler Explícito", euler),
    ("Euler Implícito", euler_implicito),
    ("RK2 Ponto Médio", rk2_ponto_medio),
    ("RK2 Heun", rk2_heun),
    ("Trapézio Implícito", trapezio_implicito),
    ("BDF2", bdf2),
    ("Adams-Bashforth 2", adams_bashforth_2),
    ("Adams-Moulton 2", adams_moulton_2),
    ("Adams-Bashforth 3", adams_bashforth_3),
    ("RK3", rk3),
    ("Adams-Moulton 3", adams_moulton_3),
    ("Adams-Bashforth 4", adams_bashforth_4),
    ("Adams-Bashforth-Moulton PC4", adams_bashforth_moulton_pc4),
    ("Dormand-Prince 5(4)", dormand_prince_54),
    ("RKF45", rkf45),
]


for i, edo in enumerate(edos, start=1):

    titulo = edo["titulo"]
    print(f"--- EDO {titulo} ---")

    expr_edo = edo["expr_edo"]
    condicoes = edo["condicoes_iniciais"]
    tf = edo["x_val"]
    h = edo["passo"]
    n = int(tf / h)

    f = edo_para_funcao(expr_edo, y(x), x)

    y0 = float(list(condicoes.values())[0])
    t0 = 0

    for nome_metodo, metodo_func in metodos:
        print(f"> Método: {nome_metodo}")

        params = inspect.signature(metodo_func).parameters

        nome_arquivo_grafico = f"{limpar_nome(titulo)}_{limpar_nome(nome_metodo)}.png"
        nome_arquivo_log = f"{limpar_nome(titulo)}_{limpar_nome(nome_metodo)}.txt"

        with open(nome_arquivo_log, "w") as f_log, redirect_stdout(f_log):

            if metodo_func == bdf2:
                val_f = f(t0, y0)
                if isinstance(val_f, np.ndarray):
                    val_f = val_f.item()
                y0_float = float(y0)
                y1_float = float(y0_float + h * val_f)
                resultados = medir_funcao(metodo_func, f, y0_float, y1_float, t0, h, n)
            elif 'n' in params:
                resultados = medir_funcao(metodo_func, f, float(y0), t0, h, n)
            elif 'tf' in params:
                resultados = medir_funcao(metodo_func, f, float(y0), t0, tf, h)
            else:
                raise Exception("Método desconhecido: assinatura incompatível")

        t_vals, y_vals = resultados[:2]
        x_exato, y_exato = resolver_edo(expr_edo, condicoes, 0, n, h)

        plt.figure(figsize=(12, 8))
        plt.plot(x_exato, y_exato, label='Solução Exata', color='blue', linestyle='--', marker='o', linewidth=1)
        plt.plot(t_vals, y_vals, label=nome_metodo, color='red', linestyle='--', marker='x', markersize=1)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title(titulo)
        plt.legend()
        plt.grid(True)

        plt.savefig(nome_arquivo_grafico)
        plt.clf()