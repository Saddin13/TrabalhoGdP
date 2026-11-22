import numpy as np

from nyquist import lugar_das_raizes
from nyquist.sistemas import SISTEMAS


def test_ramos_partem_dos_polos_de_malha_aberta():
    g = SISTEMAS["A"]
    k, ramos = lugar_das_raizes(g, k_max=100)
    assert np.allclose(np.sort_complex(ramos[0]), np.sort_complex(g.polos()))


def test_ramos_continuos():
    k, ramos = lugar_das_raizes(SISTEMAS["A"], k_max=10)
    saltos = np.max(np.abs(np.diff(ramos, axis=0)))
    assert saltos < 0.5


def test_cruza_eixo_no_ganho_critico():
    # K multiplica o numerador 10: o ganho total crítico 60 corresponde a K = 6
    k, ramos = lugar_das_raizes(SISTEMAS["A"], k_max=10)
    instavel = np.any(ramos.real > 0, axis=1)
    k_cruz = k[np.argmax(instavel)]
    assert abs(k_cruz - 6) < 0.05
