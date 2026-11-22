"""Margens de ganho e de fase a partir da resposta em frequência."""
from dataclasses import dataclass
import numpy as np


@dataclass
class Margens:
    mg: float          # margem de ganho (linear); inf se não há cruzamento de fase
    w_pc: float        # frequência de cruzamento de fase (rad/s)
    mf: float          # margem de fase (graus); inf se não há cruzamento de ganho
    w_gc: float        # frequência de cruzamento de ganho (rad/s)

    @property
    def mg_db(self):
        return 20 * np.log10(self.mg) if np.isfinite(self.mg) else np.inf

    @property
    def atraso_critico(self):
        """τcrit = MF (rad) / ωgc: maior atraso adicional tolerado."""
        if not np.isfinite(self.mf) or self.mf <= 0:
            return np.nan if np.isfinite(self.mf) else np.inf
        return np.radians(self.mf) / self.w_gc


def _bissecao(f, a, b, iters=80):
    fa = f(a)
    for _ in range(iters):
        m = np.sqrt(a * b)  # bisseção em escala logarítmica
        fm = f(m)
        if np.sign(fm) == np.sign(fa):
            a, fa = m, fm
        else:
            b = m
    return np.sqrt(a * b)


def margens(g, w_min=1e-4, w_max=1e4, n=20000):
    w = np.logspace(np.log10(w_min), np.log10(w_max), n)
    L = g(1j * w)

    # cruzamentos de ganho: |L| = 1
    mag = lambda x: np.abs(g(1j * x)) - 1
    m = np.abs(L) - 1
    idx = np.where(np.sign(m[:-1]) != np.sign(m[1:]))[0]
    mf, w_gc = np.inf, np.nan
    for i in idx:
        wc = _bissecao(mag, w[i], w[i + 1])
        fase = np.degrees(np.angle(g(1j * wc)))
        pm = (fase + 180 + 180) % 360 - 180
        if abs(pm) < abs(mf):
            mf, w_gc = pm, wc

    # cruzamentos de fase: Im(L) = 0 com Re(L) < 0
    im = lambda x: g(1j * x).imag
    idx = np.where((np.sign(L.imag[:-1]) != np.sign(L.imag[1:])))[0]
    mg, w_pc = np.inf, np.nan
    # ω = 0 também é cruzamento quando L(0) é real negativo (ex.: planta instável),
    # caso em que a margem indica a redução de ganho tolerada (MG < 1)
    if len(g.polos_no_eixo()) == 0 or np.all(np.abs(g.polos_no_eixo()) > 0):
        L0 = g(0j)
        if np.isfinite(L0) and L0.real < 0 and abs(L0.imag) < 1e-12:
            mg, w_pc = 1 / abs(L0), 0.0
    for i in idx:
        wc = _bissecao(im, w[i], w[i + 1])
        val = g(1j * wc)
        if val.real < 0:
            k = 1 / abs(val)
            if abs(np.log(k)) < abs(np.log(mg)):
                mg, w_pc = k, wc
    return Margens(mg, w_pc, mf, w_gc)
