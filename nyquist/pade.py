"""Aproximação racional de Padé para e^(−τs)."""
from math import factorial
import numpy as np


def pade(tau, ordem=1):
    """Retorna (num, den) da aproximação de Padé de ordem `ordem`, em ordem decrescente."""
    n = ordem
    c = [factorial(2 * n - k) * factorial(n) / (factorial(2 * n) * factorial(k) * factorial(n - k))
         for k in range(n + 1)]
    den = np.array([c[k] * tau ** k for k in range(n + 1)])[::-1]
    num = np.array([c[k] * (-tau) ** k for k in range(n + 1)])[::-1]
    return num, den
