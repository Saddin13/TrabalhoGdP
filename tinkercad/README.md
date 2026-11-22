# Simulação no Tinkercad

Coloque aqui os arquivos do circuito:

- `Gerencia de projetos.brd` (exportação do Tinkercad)
- `schematic.pdf` / `schematic.png` (vista esquemática)
- `desligado.png` (circuito montado)
- prints do Serial Monitor: `n=0-15s.png`, `n=10-15s.png` … `n=16-15s.png`, `n=44-15s.png` …

## Circuito

| Componente | Valor | Ligação |
|---|---|---|
| R1 | 10 kΩ | D9 → nó V1 |
| C1 | 10 µF | (+) em V1, (−) no GND |
| R2 | 100 kΩ | V1 → nó V2 |
| C2 | 1 µF | (+) em V2, (−) no GND |
| Fio | | V2 → A0 |

Planta: G(s) = 100 / (s² + 21s + 100). Firmware em `../firmware/pi_atraso/pi_atraso.ino`.

## Resultados

| N | τ (ms) | Comportamento |
|---|---|---|
| 0 | 0 | estável em 2,49 V |
| 10, 12, 13, 14 | 100–140 | estável (2,48–2,50 V) |
| **15** | **150** | **oscilação sustentada (2,19–2,78 V)** |
| 16 | 160 | oscilação maior (1,27–3,73 V) |
| 44–50 | 440–500 | oscilação saturada (0,1–4,8 V) |

Previsto por Nyquist: τcrit = MF/ωgc = 89,6° / 9,53 rad/s ≈ 164 ms (N ≈ 16).
Observado: 150 ms (N = 15), 8,7% abaixo do previsto.
