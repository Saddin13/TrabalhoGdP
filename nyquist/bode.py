"""Diagrama de Bode: módulo (dB) e fase (graus) em função de ω."""
import numpy as np


def bode(g, w_min=1e-3, w_max=1e3, n=2000):
    w = np.logspace(np.log10(w_min), np.log10(w_max), n)
    L = g(1j * w)
    mag_db = 20 * np.log10(np.abs(L))
    fase = np.degrees(np.unwrap(np.angle(L)))
    return w, mag_db, fase
