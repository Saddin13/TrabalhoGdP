"""Contorno de Nyquist e contagem de envoltórios (winding number)."""
import numpy as np

LIMIAR_FASE = np.pi / 8


def _segmentos(g, w_min=1e-4, w_max=1e4, eps=1e-6):
    """Divide o eixo imaginário em trechos retos e semicircunferências de raio eps.

    Cada segmento é (tipo, a, b): 'reta' com ω de a até b, ou 'arco' em torno
    do polo jω0 com θ de −π/2 até π/2 (desvio pela direita).
    """
    p_eixo = sorted({round(p.imag, 12) for p in g.polos_no_eixo()})
    escala = max([1.0] + [abs(p) for p in g.polos()] + [abs(z) for z in g.zeros()])
    w_max = max(w_max, 100 * escala)
    cortes = [-w_max] + p_eixo + [w_max]
    segs = []
    for i in range(len(cortes) - 1):
        a, b = cortes[i], cortes[i + 1]
        if i > 0:
            segs.append(("arco", cortes[i], None))
            a = a + eps
        if i < len(cortes) - 2:
            b = b - eps
        segs.append(("reta", a, b))
    return segs, w_min, eps


def _pontos_reta(a, b, w_min, n):
    """Pontos ω em [a, b], densos em escala logarítmica perto de 0."""
    if a < 0 < b:
        neg = -np.logspace(np.log10(-a), np.log10(w_min), n)
        pos = np.logspace(np.log10(w_min), np.log10(b), n)
        return np.concatenate([neg, [0.0], pos])
    if b <= 0 and a < 0:
        lo, hi = max(-b, w_min * 1e-3), -a
        return -np.logspace(np.log10(hi), np.log10(lo), n)
    lo, hi = max(a, w_min * 1e-3), b
    return np.logspace(np.log10(lo), np.log10(hi), n)


def contorno_nyquist(g, n=2000, refinar=True, max_iter=30):
    """Retorna (s, L) ao longo do contorno de Nyquist, percorrido no sentido horário.

    O contorno sobe pelo eixo imaginário de −j∞ a +j∞, desviando dos polos sobre
    o eixo por semicircunferências pela direita, e fecha pela semicircunferência
    de raio infinito, onde L(s) é constante para G própria (não é amostrada).
    """
    segs, w_min, eps = _segmentos(g)
    blocos = []
    for tipo, a, b in segs:
        if tipo == "reta":
            t = _pontos_reta(a, b, w_min, n)
            s = 1j * t
        else:
            t = np.linspace(-np.pi / 2, np.pi / 2, 200)
            s = 1j * a + eps * np.exp(1j * t)
        blocos.append([tipo, a, t, s])

    for bloco in blocos:
        tipo, a, t, s = bloco
        if not refinar:
            continue
        for _ in range(max_iter):
            F = 1 + g(s)
            dfase = np.abs(np.diff(np.unwrap(np.angle(F))))
            ruins = np.where(dfase > LIMIAR_FASE)[0]
            if len(ruins) == 0:
                break
            meio = (t[ruins] + t[ruins + 1]) / 2
            t = np.sort(np.concatenate([t, meio]))
            s = 1j * t if tipo == "reta" else 1j * a + eps * np.exp(1j * t)
        bloco[2], bloco[3] = t, s

    s = np.concatenate([b[3] for b in blocos])
    return s, g(s)


def numero_de_voltas(L, ponto=-1.0):
    """N: envoltórios horários de `ponto` pela curva L (anti-horários contam negativo)."""
    F = np.asarray(L) - ponto
    dtheta = np.unwrap(np.angle(F))
    return int(np.rint(-(dtheta[-1] - dtheta[0]) / (2 * np.pi)))


def criterio_nyquist(g, **kw):
    """Aplica Z = N + P. Retorna dicionário com N, P, Z e o veredito."""
    s, L = contorno_nyquist(g, **kw)
    N = numero_de_voltas(L)
    P = g.polos_instaveis()
    Z = N + P
    return {"N": N, "P": P, "Z": Z, "estavel": Z == 0, "s": s, "L": L}
