"""Gera os resultados das Tabelas 1 e 2 e da Seção 7 do artigo."""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from nyquist import (FuncaoTransferencia, criterio_nyquist, margens, tabela_routh,  # noqa: E402
                     trocas_de_sinal, ganho_critico_routh, ganho_critico_lugar, pade, raizes)
from nyquist.sistemas import SISTEMAS, D, MALHA_TINKERCAD, PLANTA_RC  # noqa: E402

SAIDA = Path(__file__).resolve().parents[1] / "resultados"


def fmt(x, nd=3):
    return "inf" if not np.isfinite(x) else round(float(x), nd)


def tabela1():
    linhas = {}
    for nome, g in SISTEMAS.items():
        r = criterio_nyquist(g)
        m = margens(g)
        linha = {"P": r["P"], "N": r["N"], "Z": r["Z"], "estavel_nyquist": bool(r["estavel"]),
                 "MG": fmt(m.mg), "MG_dB": fmt(m.mg_db, 1), "w_pc": fmt(m.w_pc),
                 "MF": fmt(m.mf, 1), "w_gc": fmt(m.w_gc)}
        if g.atraso == 0:
            cl = g.polinomio_malha_fechada()
            tab = tabela_routh(cl)
            linha["polinomio_mf"] = cl.tolist()
            linha["routh_primeira_coluna"] = [round(l[0], 3) for l in tab]
            linha["trocas_routh"] = trocas_de_sinal(tab)
            linha["polos_mf"] = [str(np.round(p, 2)) for p in raizes(cl)]
        else:
            num, den = pade(g.atraso, 1)
            gp = FuncaoTransferencia(np.polymul(g.num, num), np.polymul(g.den, den))
            linha["trocas_routh_pade1"] = trocas_de_sinal(tabela_routh(gp.polinomio_malha_fechada()))
            sem = FuncaoTransferencia(g.num, g.den)
            m0 = margens(sem)
            linha["MF_sem_atraso"] = fmt(m0.mf, 2)
            linha["w_gc_sem_atraso"] = fmt(m0.w_gc)
            linha["tau_crit_sem_atraso"] = fmt(m0.atraso_critico)
        linhas[nome] = linha
    return linhas


def tabela2():
    """Ganho crítico do Sistema D: Nyquist exato x Routh com Padé."""
    base = FuncaoTransferencia(D.num / D.num[-1], D.den, D.atraso)  # ganho unitário
    k_exato = margens(base).mg
    out = {"Nyquist (exato)": {"ordem": None, "Kcrit": fmt(k_exato, 4), "erro_%": 0.0}}
    for ordem in (1, 2):
        num, den = pade(D.atraso, ordem)
        gp = FuncaoTransferencia(np.polymul(base.num, num), np.polymul(base.den, den))
        k = ganho_critico_routh(gp)
        out[f"Routh + Padé {ordem}"] = {"ordem": len(gp.den) - 1, "Kcrit": fmt(k, 4),
                                        "erro_%": fmt(100 * (k - k_exato) / k_exato, 2)}
        k_rl = ganho_critico_lugar(gp)
        out[f"Lugar das Raízes + Padé {ordem}"] = {"ordem": len(gp.den) - 1, "Kcrit": fmt(k_rl, 4),
                                                   "erro_%": fmt(100 * (k_rl - k_exato) / k_exato, 2)}
    return out


def tinkercad():
    m = margens(MALHA_TINKERCAD)
    r = criterio_nyquist(MALHA_TINKERCAD)
    return {"planta_den": PLANTA_RC.den.tolist(), "polos_planta": [str(np.round(p, 2)) for p in PLANTA_RC.polos()],
            "N": r["N"], "P": r["P"], "Z": r["Z"],
            "MF": fmt(m.mf, 2), "w_gc": fmt(m.w_gc), "MG": fmt(m.mg),
            "tau_crit_s": fmt(m.atraso_critico, 4), "N_crit_T10ms": fmt(m.atraso_critico / 0.01, 1),
            "tau_crit_observado_s": 0.150,
            "diferenca_%": fmt(100 * (0.150 - m.atraso_critico) / m.atraso_critico, 1)}


if __name__ == "__main__":
    SAIDA.mkdir(exist_ok=True)
    res = {"tabela1": tabela1(), "tabela2": tabela2(), "tinkercad": tinkercad()}
    (SAIDA / "resultados.json").write_text(json.dumps(res, indent=2, ensure_ascii=False))
    print(json.dumps(res, indent=2, ensure_ascii=False))
