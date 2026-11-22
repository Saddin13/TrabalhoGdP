import control as ct
import numpy as np
import pytest

from nyquist import FuncaoTransferencia, margens
from nyquist.sistemas import SISTEMAS, MALHA_TINKERCAD


@pytest.mark.parametrize("nome", ["A", "B", "C"])
def test_margens_iguais_ao_python_control(nome):
    g = SISTEMAS[nome]
    gm, pm, wpc, wgc = ct.margin(ct.tf(g.num, g.den))
    m = margens(g)
    if np.isfinite(gm):
        assert m.mg == pytest.approx(gm, rel=1e-4)
        assert m.w_pc == pytest.approx(wpc, rel=1e-4)
    else:
        assert not np.isfinite(m.mg)
    if np.isfinite(pm):
        assert m.mf == pytest.approx(pm, rel=1e-4)
        assert m.w_gc == pytest.approx(wgc, rel=1e-4)
    else:
        assert not np.isfinite(m.mf)


def test_margens_com_atraso():
    g = SISTEMAS["D"]
    m = margens(g)
    # com atraso, fase = fase sem atraso − ωτ
    sem = FuncaoTransferencia(g.num, g.den)
    m0 = margens(sem)
    assert m.w_gc == pytest.approx(m0.w_gc, rel=1e-6)
    assert m.mf == pytest.approx(m0.mf - np.degrees(0.5 * m0.w_gc), rel=1e-4)


def test_atraso_critico_tinkercad():
    m = margens(MALHA_TINKERCAD)
    assert m.atraso_critico == pytest.approx(0.164, abs=1e-3)
