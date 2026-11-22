import numpy as np
import pytest

from nyquist import (FuncaoTransferencia, tabela_routh, trocas_de_sinal, estavel_routh,
                     ganho_critico_routh, ganho_critico_lugar, margens, pade)
from nyquist.sistemas import SISTEMAS


@pytest.mark.parametrize("polos", [
    [-1, -2, -3], [1, -2, -3], [1, 2, -3], [-1 + 2j, -1 - 2j, -4],
    [0.5 + 1j, 0.5 - 1j, -1, -2, -7], [-1, -2, -3, -4, -5, -6],
])
def test_trocas_de_sinal_igual_raizes_instaveis(polos):
    c = np.poly(polos).real
    assert trocas_de_sinal(tabela_routh(c)) == int(np.sum(np.array(polos).real > 0))


def test_zero_na_primeira_coluna():
    # s^4 + s^3 + 2s^2 + 2s + 3: zero na 1ª coluna, 2 raízes instáveis
    c = [1, 1, 2, 2, 3]
    assert trocas_de_sinal(tabela_routh(c)) == int(np.sum(np.roots(c).real > 0))


def test_linha_nula():
    # s^3 + s^2 + 4s + 4 = (s+1)(s^2+4): raízes ±2j, linha nula
    tab = tabela_routh([1, 1, 4, 4])
    assert trocas_de_sinal(tab) == 0


def test_sistema_a_routh():
    tab = tabela_routh(SISTEMAS["A"].polinomio_malha_fechada())
    assert [l[0] for l in tab] == pytest.approx([1, 6, 50 / 6, 16])


@pytest.mark.parametrize("nome,kcrit", [("A", 6.0), ("C", 60 / 7 / 5)])
def test_ganho_critico_routh_e_lugar(nome, kcrit):
    g = SISTEMAS[nome]
    base = FuncaoTransferencia(g.num, g.den)
    assert ganho_critico_routh(base) == pytest.approx(kcrit, rel=1e-6)
    assert ganho_critico_lugar(base) == pytest.approx(kcrit, rel=1e-6)
    assert margens(base).mg == pytest.approx(kcrit, rel=1e-4)


def test_pade_coeficientes():
    num, den = pade(0.5, 1)
    assert np.allclose(num, [-0.25, 1]) and np.allclose(den, [0.25, 1])
    num, den = pade(0.5, 2)
    assert np.allclose(den, [0.25 / 12, 0.25, 1]) and np.allclose(num, [0.25 / 12, -0.25, 1])


def test_pade_igual_python_control():
    import control as ct
    for ordem in (1, 2, 3):
        n_ref, d_ref = ct.pade(0.5, ordem)
        num, den = pade(0.5, ordem)
        assert np.allclose(num / den[0], np.array(n_ref) / d_ref[0])
        assert np.allclose(den / den[0], np.array(d_ref) / d_ref[0])
