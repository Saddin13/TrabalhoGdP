"""Representação de funções de transferência G(s) = N(s)/D(s)·e^(−τs)."""
import numpy as np


def horner(coef, s):
    """Avalia um polinômio (coeficientes em ordem decrescente) pelo esquema de Horner.

    Usa n multiplicações e n somas para grau n. Aceita s escalar ou array.
    """
    s = np.asarray(s, dtype=complex)
    y = np.zeros_like(s)
    for c in coef:
        y = y * s + c
    return y


def raizes(coef):
    """Raízes de um polinômio como autovalores da matriz companheira."""
    c = np.trim_zeros(np.asarray(coef, dtype=float), "f")
    if len(c) <= 1:
        return np.array([], dtype=complex)
    c = c / c[0]
    n = len(c) - 1
    comp = np.zeros((n, n))
    comp[0, :] = -c[1:]
    comp[1:, :-1] = np.eye(n - 1)
    return np.linalg.eigvals(comp).astype(complex)


class FuncaoTransferencia:
    """G(s) = N(s)/D(s)·e^(−τs), coeficientes em ordem decrescente de potência."""

    def __init__(self, num, den, atraso=0.0):
        self.num = np.trim_zeros(np.atleast_1d(np.asarray(num, dtype=float)), "f")
        self.den = np.trim_zeros(np.atleast_1d(np.asarray(den, dtype=float)), "f")
        if len(self.den) == 0:
            raise ValueError("denominador nulo")
        if len(self.num) > len(self.den):
            raise ValueError("função de transferência imprópria")
        self.atraso = float(atraso)

    def avaliar(self, s):
        s = np.asarray(s, dtype=complex)
        g = horner(self.num, s) / horner(self.den, s)
        if self.atraso:
            g = g * np.exp(-self.atraso * s)
        return g

    __call__ = avaliar

    def polos(self):
        return raizes(self.den)

    def zeros(self):
        return raizes(self.num)

    def polos_instaveis(self, tol=1e-9):
        """P: número de polos no semiplano direito aberto."""
        return int(np.sum(self.polos().real > tol))

    def polos_no_eixo(self, tol=1e-9):
        p = self.polos()
        return p[np.abs(p.real) <= tol]

    def serie(self, outra):
        """Produto G1(s)·G2(s)."""
        return FuncaoTransferencia(np.polymul(self.num, outra.num),
                                   np.polymul(self.den, outra.den),
                                   self.atraso + outra.atraso)

    def com_ganho(self, k):
        return FuncaoTransferencia(k * self.num, self.den, self.atraso)

    def polinomio_malha_fechada(self, k=1.0):
        """D(s) + K·N(s) (realimentação unitária negativa, sem atraso)."""
        if self.atraso:
            raise ValueError("polinômio característico não existe com atraso: use Padé")
        return np.polyadd(self.den, k * self.num)

    def __repr__(self):
        a = f", atraso={self.atraso}" if self.atraso else ""
        return f"FuncaoTransferencia(num={self.num.tolist()}, den={self.den.tolist()}{a})"
