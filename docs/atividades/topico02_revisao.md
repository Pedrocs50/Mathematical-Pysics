# Tópico 2: revisão item por item

[Índice](../README.md) · [Como explicar](topico02_como_explicar.md) · [Código](../../topic02_variational/)

Revisão de 25/11/2026 contra o enunciado (apostila, páginas 32 a 34). No mesmo dia, os problemas encontrados foram corrigidos, os scripts rodados de novo e os notebooks salvos com as saídas.

Legenda: ✅ atende

## O que foi corrigido

| Item | Problema encontrado | Correção |
|---|---|---|
| 3.1 | faltava o desenho do pêndulo | desenho em matplotlib, salvo em `outputs/h3_desenho.png` |
| 2.5 | ẋ(t) não era plotado com atrito | a figura agora tem x, ẋ e energia |
| 2.5 | a figura `h2_com_atrito.png` não correspondia ao script | gerada de novo |
| 2.5 | o regime de amortecimento só valia para o poço simples | tabela de γ / (2√(mV'')) para os dois poços |
| 2.6 | o cruzamento com atrito era calculado e não aparecia | novo mapa, contagem e discussão |
| 2.6 | partir de x₀ = 0 contava como cruzamento | agora exige troca estrita de sinal |
| 3.6 | "em 90° o período real é ~15% maior que T₀" | corrigido para ~18% (15,3% é quanto T₀ subestima T) |
| 1.6 | dedução curta | definição de δS como d/dε e termo de borda explícito |
| notebooks | sem saídas; as figuras não apareceriam na célula | backend ajustado, notebooks executados e salvos |

## Situação item por item

| Item | Situação | Observação |
|---|---|---|
| 1.1 a 1.5 | ✅ | a rotina é conferida contra `euler_equations` |
| 1.6 | ✅ | reescrever com as próprias palavras |
| 2.1 a 2.4 | ✅ | desvio de energia cerca de 4 × 10⁻¹⁰ |
| 2.5 | ✅ | |
| 2.6 | ✅ | com γ = 0,2, 1218 de 1420 condições com E₀ > 0,25 cruzam |
| 2.7 | ✅ | reescrever; dá para enriquecer com o [contraste entre máximo e mínimo de V](../teoria/02_calculo_variacional.md#8-estacionário-ou-mínimo) |
| 3.1 a 3.3 | ✅ | |
| 3.4 | ✅ | dedução breve; dá para mostrar a solução θ₀ cos ω₀t |
| 3.5 e 3.6 | ✅ | |
| 3.7 | ✅ | reescrever; dá para mostrar dH/dt = 0 |
| 3.8 | ✅ | erro de 4 × 10⁻⁴ no balanço, por diferenças finitas |

## Antes de entregar

- [ ] Reescrever 1.6, 2.7 e 3.7 com as próprias palavras.
- [ ] Ler as saídas dos três notebooks e saber explicar cada figura.
- [ ] Confirmar com o professor o formato e o prazo de entrega.

> Gerar um notebook com `tools/py_to_ipynb.py` apaga as saídas dele: depois de gerar, execute o notebook e salve. **No Tópico 1, não gere:** as respostas só existem no notebook.
