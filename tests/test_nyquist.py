import control as ct
import numpy as np
import pytest

from nyquist import FuncaoTransferencia, criterio_nyquist
from nyquist.sistemas import SISTEMAS, MALHA_TINKERCAD


def para_control(g):
    sys_ = ct.tf(g.num, g.den)
    return sys_


@pytest.mark.parametrize("nome", ["A", "B", "C"])
def test_z_igual_polos_instaveis_de_malha_fechada(nome):
    g = SISTEMAS[nome]
    r = criterio_nyquist(g)
    mf = ct.feedback(para_control(g), 1)
    z_ref = int(np.sum(ct.poles(mf).real > 0))
    assert r["Z"] == z_ref


@pytest.mark.parametrize("nome", ["A", "B", "C"])
def test_envoltorios_iguais_ao_python_control(nome):
    g = SISTEMAS[nome]
    resp = ct.nyquist_response(para_control(g))
    # mesma convenção: N horário, Z = N + P
    assert criterio_nyquist(g)["N"] == int(resp.count)


def test_sistema_b_envolve_uma_vez_anti_horario():
    r = criterio_nyquist(SISTEMAS["B"])
    assert (r["P"], r["N"], r["Z"]) == (1, -1, 0)


@pytest.mark.parametrize("k,esperado", [(10, 0), (59, 0), (61, 2), (200, 2)])
def test_sistema_a_instavel_acima_do_ganho_critico(k, esperado):
    g = FuncaoTransferencia([k], np.poly([-1, -2, -3]))
    assert criterio_nyquist(g)["Z"] == esperado


def test_integrador_no_eixo_imaginario():
    # PI·planta RC tem polo em s = 0: o contorno precisa desviar dele
    r = criterio_nyquist(MALHA_TINKERCAD)
    assert r["P"] == 0 and r["Z"] == 0


def test_polos_imaginarios_puros():
    # G = 1/(s(s^2+1)) com K pequeno: malha fechada s^3+s+K tem 2 polos instáveis
    g = FuncaoTransferencia([0.5], [1, 0, 1, 0])
    z_ref = int(np.sum(np.roots([1, 0, 1, 0.5]).real > 0))
    assert criterio_nyquist(g)["Z"] == z_ref


@pytest.mark.parametrize("tau,estavel", [(0.5, True), (1.0, True), (1.5, False), (3.0, False)])
def test_atraso_exato(tau, estavel):
    # Sistema D sem atraso tem τcrit = 1,22 s
    g = FuncaoTransferencia([2], np.polymul([1, 1], [0.5, 1]), atraso=tau)
    assert criterio_nyquist(g)["estavel"] == estavel
