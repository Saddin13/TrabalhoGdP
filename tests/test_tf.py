import numpy as np
import pytest

from nyquist import FuncaoTransferencia, horner, raizes


def test_horner_igual_polyval():
    c = [3.0, -2.0, 0.5, 7.0]
    s = np.array([0.3 + 1j, -2.0, 5j])
    assert np.allclose(horner(c, s), np.polyval(c, s))


def test_raizes_companheira():
    c = np.poly([-1, -2, -3, 4 + 1j, 4 - 1j]).real
    assert np.allclose(np.sort_complex(raizes(c)), np.sort_complex(np.roots(c)))


def test_atraso_nao_altera_modulo():
    g0 = FuncaoTransferencia([1], [1, 1])
    g = FuncaoTransferencia([1], [1, 1], atraso=0.7)
    w = np.logspace(-2, 2, 50)
    assert np.allclose(np.abs(g(1j * w)), np.abs(g0(1j * w)))
    assert np.allclose(np.angle(g(1j * w) / g0(1j * w)), np.angle(np.exp(-0.7j * w)))


def test_impropria_rejeitada():
    with pytest.raises(ValueError):
        FuncaoTransferencia([1, 0, 0], [1, 1])
