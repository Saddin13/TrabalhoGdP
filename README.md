# Análise de Estabilidade pelo Critério de Nyquist

Ferramenta em Python para analisar a estabilidade de sistemas lineares invariantes no tempo
pelo Critério de Nyquist, comparada com Routh-Hurwitz, Lugar das Raízes e Diagrama de Bode.
O núcleo (`nyquist/`) usa apenas NumPy; a biblioteca `python-control` é usada só como
referência nos testes e no benchmark.

## Estrutura

```
nyquist/            núcleo da ferramenta
  tf.py             FuncaoTransferencia: N(s)/D(s)·e^(−τs), Horner, raízes por matriz companheira
  nyquist.py        contorno de Nyquist (com desvio de polos no eixo jω), winding number, Z = N + P
  margens.py        margens de ganho e de fase, atraso crítico τcrit = MF/ωgc
  routh.py          tabela de Routh (zero na 1ª coluna e linha nula), ganho crítico
  lugar_raizes.py   Lugar das Raízes com passo adaptativo e ganho crítico
  bode.py           módulo (dB) e fase
  pade.py           aproximação de Padé de qualquer ordem
  sistemas.py       Sistemas A, B, C, D e a malha do Tinkercad
experimentos/
  resultados.py     Tabelas 1 e 2 e resultados da Seção 7   → resultados/resultados.json
  escalabilidade.py Tabela 3 (tempo de CPU e erro, n = 5, 10, 20) → resultados/escalabilidade.json
  figuras.py        Figuras 1, 2 e 4 → figuras/
tests/              testes unitários (validação cruzada com python-control)
firmware/           PI com atraso simulado para Arduino (Tinkercad)
tinkercad/          circuito, esquemático e prints da simulação
```

## Como reproduzir

```bash
python -m venv .venv && source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r requirements.txt
python -m pytest -q
python experimentos/resultados.py
python experimentos/escalabilidade.py
python experimentos/figuras.py
```

## Exemplo

```python
from nyquist import FuncaoTransferencia, criterio_nyquist, margens

g = FuncaoTransferencia([2], [0.5, 1.5, 1], atraso=0.5)   # 2e^(−0,5s)/((s+1)(0,5s+1))
r = criterio_nyquist(g)
print(r["N"], r["P"], r["Z"], r["estavel"])               # 0 0 0 True
m = margens(g)
print(m.mg, m.mf, m.atraso_critico)
```

## Convenções

- Coeficientes de polinômios em ordem decrescente de potência.
- N conta envoltórios **horários** de (−1, 0); anti-horários contam como negativos. Estável ⇔ Z = N + P = 0.
- Realimentação unitária negativa.
