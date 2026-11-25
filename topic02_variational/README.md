# Tópico 2: cálculo variacional (Prof. Reinaldo)

As três atividades estão na seção "Hands-on Activities" (p. 32-34) da apostila `variational_principle.pdf`, enviada por email em 16/11. Os itens abaixo são um resumo; vale o texto da apostila.

| Atividade | Arquivos | Figuras em `outputs/` |
|---|---|---|
| 1. Euler-Lagrange simbólico | `handson_1_euler_lagrange.py` / `.ipynb` | – |
| 2. Poços simples e duplo | `handson_2_potential_wells.py` / `.ipynb` | `h2_*.png` |
| 3. Pêndulo simples | `handson_3_pendulum.py` / `.ipynb` | `h3_*.png` |

Explicação dos resultados: [como explicar](../docs/atividades/topico02_como_explicar.md). O que foi corrigido e o que fazer antes de entregar: [revisão](../docs/atividades/topico02_revisao.md).

Código usado: `cap422/variational.py` e `cap422/dynamics.py`. Testes: `tests/test_variational.py`.

## Hands-On 1: Euler-Lagrange simbólico com SymPy

- **1.1** Definir x, y(x) e y'(x) simbolicamente.
- **1.2** Construir a ação S[y] e dizer o que é função e o que é funcional.
- **1.3** Implementar uma rotina que calcula d/dx(∂L/∂y') − ∂L/∂y para uma L qualquer.
- **1.4** Testar com L = ½ y'² + ½ ω² y², mostrar a equação e pedir a solução simbólica.
- **1.5** Repetir com L = ½ m y'² − ½ k y² e comparar os comportamentos.
- **1.6** Explicar a relação entre Euler-Lagrange e a primeira variação da ação.

## Hands-On 2: movimento em poços simples e duplo

- **2.1** Implementar V1 = ½ k x² e V2 = a (x² − b²)², plotar e achar os pontos estacionários.
- **2.2** Derivar V'(x) e a equação de Euler-Lagrange de cada um com SymPy.
- **2.3** Escrever como sistema de primeira ordem e integrar com `scipy.integrate.solve_ivp`.
- **2.4** Plotar x(t), velocidade e energia mecânica com atrito zero e avaliar a conservação.
- **2.5** Incluir atrito de Rayleigh (R = ½ γ ẋ²) com γ = 0,2, 0,6 e 2,0.
- **2.6** No poço duplo, variar posição e velocidade iniciais, ver quando cruza a barreira e para qual mínimo vai com atrito.
- **2.7** Discutir por que um mínimo local de V não significa que a ação é mínima.

## Hands-On 3: pêndulo simples

- **3.1** Desenhar o pêndulo e identificar θ(t), m, l e g.
- **3.2** Obter T, V (zero embaixo) e L, e escrever a ação.
- **3.3** Obter a equação de Euler-Lagrange com SymPy, sem partir da força tangencial.
- **3.4** Aproximação de pequenos ângulos, frequência e período.
- **3.5** Resolver numericamente a equação não linear e a linearizada para 5, 30 e 90 graus, partindo do repouso.
- **3.6** Comparar trajetórias e períodos, quantificar o erro e verificar a conservação da energia.
- **3.7** Dizer qual simetria, pelo teorema de Noether, leva à conservação da energia.
- **3.8** Opcional: atrito de Rayleigh e balanço de energia.

## Como rodar

A partir da raiz do repositório, por exemplo:

```powershell
.\.venv\Scripts\python.exe topic02_variational/handson_1_euler_lagrange.py
```

Para regenerar um notebook a partir do script (o notebook gerado vem sem saídas; execute-o e salve):

```powershell
.\.venv\Scripts\python.exe tools/py_to_ipynb.py topic02_variational/handson_1_euler_lagrange.py
```
