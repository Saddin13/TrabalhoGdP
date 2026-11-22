"""Gera as Figuras 1, 2 e 4 do artigo em figuras/."""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from nyquist import (FuncaoTransferencia, criterio_nyquist, lugar_das_raizes, bode,  # noqa: E402
                     margens)
from nyquist.sistemas import SISTEMAS, PI, PLANTA_RC  # noqa: E402

PASTA = Path(__file__).resolve().parents[1] / "figuras"
plt.rcParams.update({"font.family": "serif", "font.size": 9, "axes.grid": True,
                     "grid.alpha": 0.3, "savefig.dpi": 200})
COR = "#1f4e79"


def _nyquist(ax, g, titulo, lim=None):
    r = criterio_nyquist(g)
    L = r["L"]
    ok = np.isfinite(L) & (np.abs(L) < 1e3)
    ax.plot(L.real[ok], L.imag[ok], color=COR, lw=1)
    ax.plot(-1, 0, "r+", ms=10, mew=1.5)
    ax.set_title(f"{titulo}\nP = {r['P']}, N = {r['N']}, Z = {r['Z']}")
    ax.set_xlabel("Re")
    ax.set_ylabel("Im")
    ax.axhline(0, color="k", lw=0.5)
    ax.axvline(0, color="k", lw=0.5)
    if lim:
        ax.set_xlim(*lim[0])
        ax.set_ylim(*lim[1])


def figura1():
    fig, axs = plt.subplots(2, 2, figsize=(7, 6))
    lims = {"A": ((-1.5, 2), (-1.5, 1.5)), "B": ((-4.5, 1), (-3, 3)),
            "C": ((-1.3, 1), (-0.8, 0.8)), "D": ((-1.5, 2.2), (-1.8, 1.8))}
    nomes = {"A": "Sistema A (padrão)", "B": "Sistema B (instável em MA)",
             "C": "Sistema C (fase não-mínima)", "D": "Sistema D (atraso 0,5 s)"}
    for ax, (k, g) in zip(axs.flat, SISTEMAS.items()):
        _nyquist(ax, g, nomes[k], lims[k])
    fig.tight_layout()
    fig.savefig(PASTA / "figura1_nyquist.png")


def figura2():
    fig, axs = plt.subplots(2, 2, figsize=(7, 6))
    for col, k in enumerate(("A", "C")):
        g = SISTEMAS[k]
        base = FuncaoTransferencia(g.num / abs(g.num[-1]), g.den)
        kk, ramos = lugar_das_raizes(base, k_max=100 if k == "A" else 15, n=800)
        ax = axs[0, col]
        for j in range(ramos.shape[1]):
            ax.plot(ramos[:, j].real, ramos[:, j].imag, color=COR, lw=1)
        ax.plot(g.polos().real, g.polos().imag, "kx")
        if len(g.zeros()):
            ax.plot(g.zeros().real, g.zeros().imag, "ko", mfc="none")
        ax.axvline(0, color="r", lw=0.6)
        ax.set_title(f"Lugar das Raízes, Sistema {k}")
        ax.set_xlabel("Re(s)")
        ax.set_ylabel("Im(s)")
        w, mag, fase = bode(g, 1e-2, 1e2)
        m = margens(g)
        ax1 = axs[1, col]
        ax1.semilogx(w, mag, color=COR, lw=1, label="|L| (dB)")
        ax1.set_ylabel("Módulo (dB)")
        ax1.set_xlabel("ω (rad/s)")
        ax2 = ax1.twinx()
        ax2.semilogx(w, fase, color="#c0504d", lw=1, ls="--", label="fase")
        ax2.axhline(-180, color="#c0504d", lw=0.5)
        ax2.set_ylabel("Fase (graus)")
        ax2.grid(False)
        ax1.axhline(0, color="k", lw=0.5)
        ax1.set_title(f"Bode, Sistema {k}: MG = {m.mg:.2f}")
    fig.tight_layout()
    fig.savefig(PASTA / "figura2_lugar_bode.png")


def figura4():
    """Nyquist da malha do Tinkercad para atrasos crescentes."""
    fig, ax = plt.subplots(figsize=(4.5, 4))
    malha = PI.serie(PLANTA_RC)
    for tau, cor in ((0.0, "#1f4e79"), (0.10, "#4f81bd"), (0.164, "#c0504d"), (0.25, "#9bbb59")):
        g = FuncaoTransferencia(malha.num, malha.den, tau)
        w = np.logspace(-0.5, 2.5, 4000)
        L = g(1j * w)
        ax.plot(L.real, L.imag, color=cor, lw=1, label=f"τ = {1000 * tau:.0f} ms")
    ax.plot(-1, 0, "r+", ms=10, mew=1.5)
    ax.set_xlim(-2, 1)
    ax.set_ylim(-2.5, 0.5)
    ax.set_xlabel("Re")
    ax.set_ylabel("Im")
    ax.set_title(f"Malha PI + planta RC (τcrit = {1000 * margens(malha).atraso_critico:.0f} ms)")
    ax.legend(loc="lower right", fontsize=7)
    fig.tight_layout()
    fig.savefig(PASTA / "figura4_tinkercad_nyquist.png")


if __name__ == "__main__":
    PASTA.mkdir(exist_ok=True)
    figura1()
    figura2()
    figura4()
    print("figuras salvas em", PASTA)
