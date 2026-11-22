"""Sistemas de referência usados no artigo."""
import numpy as np
from .tf import FuncaoTransferencia

# Sistema A: estável em malha aberta, 3ª ordem
A = FuncaoTransferencia([10], np.poly([-1, -2, -3]))
# Sistema B: instável em malha aberta (polo em s = +1)
B = FuncaoTransferencia([10, 20], np.poly([1, -5]))
# Sistema C: fase não-mínima (zero em s = +1)
C = FuncaoTransferencia([-5, 5], np.poly([-1, -2, -3]))
# Sistema D: atraso de transporte de 0,5 s
D = FuncaoTransferencia([2], np.polymul([1, 1], [0.5, 1]), atraso=0.5)

SISTEMAS = {"A": A, "B": B, "C": C, "D": D}

# Planta RC de 2 estágios do Tinkercad (R1=10k, C1=10uF, R2=100k, C2=1uF)
R1, C1, R2, C2 = 10e3, 10e-6, 100e3, 1e-6
PLANTA_RC = FuncaoTransferencia(
    [1.0], [R1 * C1 * R2 * C2, R1 * C1 + R2 * C2 + R1 * C2, 1.0])
# Controlador PI do firmware: Kp = 2, Ki = 1
KP, KI = 2.0, 1.0
PI = FuncaoTransferencia([KP, KI], [1.0, 0.0])
MALHA_TINKERCAD = PI.serie(PLANTA_RC)
