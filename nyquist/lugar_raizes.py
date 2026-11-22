"""Lugar das Raízes: trajetória dos polos de malha fechada com K de 0 a K_max."""
import numpy as np
from .tf import raizes


def _associar(anteriores, novas):
    """Ordena `novas` para que cada raiz fique no ramo da raiz anterior mais próxima."""
    dist = np.abs(anteriores[:, None] - novas[None, :])
    saida = np.empty_like(anteriores)
    livres = np.ones(len(novas), dtype=bool)
    for i in np.argsort(dist.min(axis=1)):
        j = int(np.argmin(np.where(livres, dist[i], np.inf)))
        saida[i] = novas[j]
        livres[j] = False
    return saida


def lugar_das_raizes(g, k_max=100.0, n=400, fator_min=0.1):
    """Retorna (K, ramos) com ramos[i, j] = j-ésimo polo para K[i].

    O passo em K diminui (até `fator_min` do passo nominal) quando algum ramo se
    aproxima do eixo imaginário.
    """
    ks = [0.0]
    r = raizes(g.polinomio_malha_fechada(0.0))
    ramos = [r]
    passo = k_max / n
    k = 0.0
    while k < k_max:
        dist = np.min(np.abs(ramos[-1].real)) if len(ramos[-1]) else 1.0
        h = passo * min(1.0, max(fator_min, dist))
        k = min(k + h, k_max)
        novas = raizes(g.polinomio_malha_fechada(k))
        ramos.append(_associar(ramos[-1], novas))
        ks.append(k)
    return np.array(ks), np.array(ramos)


def ganho_critico_lugar(g, k_max=1e3, n=400, iters=200):
    """Primeiro K em que um ramo cruza para o semiplano direito, refinado por bisseção."""
    def instavel(k):
        return np.any(raizes(g.polinomio_malha_fechada(k)).real > 0)

    if instavel(0.0):
        return 0.0
    ks = np.concatenate([[0.0], np.logspace(-3, np.log10(k_max), n)])
    anterior = 0.0
    for k in ks[1:]:
        if instavel(k):
            lo, hi = anterior, k
            for _ in range(iters):
                if hi - lo <= 1e-15 * hi:
                    break
                m = (lo + hi) / 2
                lo, hi = (lo, m) if instavel(m) else (m, hi)
            return (lo + hi) / 2
        anterior = k
    return np.inf
