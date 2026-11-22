"""Tabela 3: tempo de CPU e erro numérico em função da ordem do sistema."""
import json
import platform
import sys
import time
from pathlib import Path

import control as ct
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from nyquist import (FuncaoTransferencia, criterio_nyquist, margens, tabela_routh,  # noqa: E402
                     trocas_de_sinal, ganho_critico_routh, lugar_das_raizes,
                     ganho_critico_lugar, bode, horner)

SAIDA = Path(__file__).resolve().parents[1] / "resultados"
REPETICOES = 100
ORDENS = (5, 10, 20)


def sistema_aleatorio(n, rng):
    """Polos reais negativos em [−10, −0,5]; ganho = metade do ganho crítico."""
    polos = -rng.uniform(0.5, 10.0, n)
    den = np.poly(polos)
    g1 = FuncaoTransferencia([1.0], den)
    gm_ref = ct.margin(ct.tf([1.0], den))[0]
    k0 = 0.5 * gm_ref
    return FuncaoTransferencia([k0], den), polos, k0, gm_ref


def cronometrar(f, rep=REPETICOES):
    f()  # aquecimento
    t0 = time.process_time()
    for _ in range(rep):
        f()
    return 1000 * (time.process_time() - t0) / rep  # ms por execução


def erro_rel(x, ref):
    return abs(x - ref) / abs(ref)


def medir(n, rng):
    g, polos, k0, gm_ref = sistema_aleatorio(n, rng)
    kcrit_ref = gm_ref  # ganho total crítico do numerador
    base = FuncaoTransferencia([1.0], g.den)
    cl = g.polinomio_malha_fechada()
    w = np.logspace(-3, 3, 2000)

    t = {
        "Nyquist": cronometrar(lambda: criterio_nyquist(g)),
        "Routh-Hurwitz": cronometrar(lambda: trocas_de_sinal(tabela_routh(cl))),
        "Lugar das Raizes": cronometrar(lambda: lugar_das_raizes(base, k_max=2 * kcrit_ref, n=400)),
        "Bode": cronometrar(lambda: bode(g, n=2000)),
    }

    # erros em relação ao python-control / à forma fatorada
    e = {
        "Nyquist": erro_rel(margens(base).mg, kcrit_ref),
        "Routh-Hurwitz": erro_rel(ganho_critico_routh(base), kcrit_ref),
        "Lugar das Raizes": erro_rel(ganho_critico_lugar(base, k_max=10 * kcrit_ref), kcrit_ref),
    }
    _, mag_db, fase = bode(g, n=2000)
    ref = ct.frequency_response(ct.tf(g.num, g.den), w)
    mag_ref = 20 * np.log10(np.abs(ref.complex))
    e["Bode"] = float(np.max(np.abs(mag_db - mag_ref) / np.maximum(1, np.abs(mag_ref))))

    # Horner x forma fatorada D(s) = Π(s − p_i), no eixo jω
    s = 1j * w
    fatorada = np.prod([s - p for p in polos], axis=0)
    e_horner = float(np.max(np.abs(horner(g.den, s) - fatorada) / np.abs(fatorada)))

    # estabilidade dos veredictos
    z_ok = criterio_nyquist(g)["Z"] == 0 and trocas_de_sinal(tabela_routh(cl)) == 0
    return t, e, e_horner, z_ok


if __name__ == "__main__":
    rng = np.random.default_rng(2026)
    tabela = {}
    for n in ORDENS:
        t, e, eh, ok = medir(n, rng)
        tabela[n] = {"tempo_ms": {k: round(v, 3) for k, v in t.items()},
                     "erro_rel": {k: float(f"{v:.2e}") for k, v in e.items()},
                     "erro_horner": float(f"{eh:.2e}"), "veredictos_estavel": bool(ok)}
        print(n, json.dumps(tabela[n], ensure_ascii=False))
    maquina = {"processador": platform.processor() or platform.machine(),
               "sistema": platform.platform(), "python": platform.python_version(),
               "numpy": np.__version__, "python-control": ct.__version__,
               "repeticoes": REPETICOES}
    try:
        for linha in open("/proc/cpuinfo"):
            if linha.startswith("model name"):
                maquina["processador"] = linha.split(":", 1)[1].strip()
                break
    except OSError:
        pass
    SAIDA.mkdir(exist_ok=True)
    (SAIDA / "escalabilidade.json").write_text(
        json.dumps({"maquina": maquina, "tabela": tabela}, indent=2, ensure_ascii=False))
    print(json.dumps(maquina, ensure_ascii=False))
