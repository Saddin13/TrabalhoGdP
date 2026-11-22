"""Critério de Routh-Hurwitz."""
import numpy as np

EPS = 1e-9


def tabela_routh(coef, eps=EPS):
    """Monta a tabela de Routh (lista de linhas) para coeficientes em ordem decrescente.

    Casos especiais: zero na primeira coluna → substituído por eps;
    linha inteira nula → derivada do polinômio auxiliar da linha anterior.
    """
    c = [float(x) for x in np.trim_zeros(np.asarray(coef, dtype=float), "f")]
    n = len(c) - 1
    largura = n // 2 + 1
    l1 = c[0::2] + [0.0] * (largura - len(c[0::2]))
    l2 = c[1::2] + [0.0] * (largura - len(c[1::2]))
    tab = [l1, l2]
    for k in range(2, n + 1):
        a, b = tab[-2], tab[-1]
        if all(abs(x) < 1e-12 for x in b):          # linha nula
            grau = n - (k - 2)                      # grau do polinômio auxiliar
            b = [a[i] * (grau - 2 * i) for i in range(largura)]
            tab[-1] = b
        if abs(b[0]) < 1e-12:                       # zero na primeira coluna
            b = list(b)
            b[0] = eps
            tab[-1] = b
        nova = [(b[0] * a[i + 1] - a[0] * b[i + 1]) / b[0] for i in range(largura - 1)] + [0.0]
        tab.append(nova)
    return tab[: n + 1]


def trocas_de_sinal(tab):
    col = [linha[0] for linha in tab]
    return int(np.sum(np.sign(col[:-1]) != np.sign(col[1:])))


def estavel_routh(coef):
    tab = tabela_routh(coef)
    return trocas_de_sinal(tab) == 0 and all(l[0] > 0 for l in tab)


def ganho_critico_routh(g, k_max=1e30, iters=200):
    """Maior K tal que D(s) + K·N(s) continua estável (bisseção sobre o critério)."""
    if not estavel_routh(g.polinomio_malha_fechada(0.0)):
        return 0.0
    lo, hi = 0.0, 1.0
    while estavel_routh(g.polinomio_malha_fechada(hi)):
        lo, hi = hi, hi * 2
        if hi > k_max:
            return np.inf
    for _ in range(iters):
        m = (lo + hi) / 2
        if hi - lo <= 1e-15 * hi:
            break
        if estavel_routh(g.polinomio_malha_fechada(m)):
            lo = m
        else:
            hi = m
    return (lo + hi) / 2
