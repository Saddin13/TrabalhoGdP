"""Ferramenta de análise de estabilidade pelo Critério de Nyquist.

Núcleo implementado apenas com NumPy, sem toolboxes de controle.
"""
from .tf import FuncaoTransferencia, horner, raizes
from .nyquist import contorno_nyquist, numero_de_voltas, criterio_nyquist
from .margens import margens, Margens
from .routh import tabela_routh, trocas_de_sinal, estavel_routh, ganho_critico_routh
from .lugar_raizes import lugar_das_raizes, ganho_critico_lugar
from .bode import bode
from .pade import pade

__all__ = [
    "FuncaoTransferencia", "horner", "raizes",
    "contorno_nyquist", "numero_de_voltas", "criterio_nyquist",
    "margens", "Margens",
    "tabela_routh", "trocas_de_sinal", "estavel_routh", "ganho_critico_routh",
    "lugar_das_raizes", "ganho_critico_lugar",
    "bode", "pade",
]
