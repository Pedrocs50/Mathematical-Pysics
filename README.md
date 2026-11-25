# CAP-422: Física Matemática Computacional (INPE)

Repositório de apoio da disciplina: Hands-On computacionais com SymPy e NumPy, rodando 100% local no VS Code (sem Colab).

## Como começar

```bash
# dentro da pasta do projeto
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

No VS Code: abra a pasta, instale as extensões recomendadas (Python e Jupyter) e escolha o interpretador `.venv` (Ctrl+Shift+P → *Python: Select Interpreter*).

## Hands-On 1: tensor de condutividade 2D

| Arquivo | O que é |
|---|---|
| `hands_on_01.ipynb` | Notebook completo, com as respostas das 5 perguntas de discussão. |
| `hands_on_01.py` | Mesma coisa em script com células `# %%` (rode célula por célula ou `python hands_on_01.py`). |
| `cap422/tensors2d.py` | Funções SymPy/NumPy reutilizáveis (σ, R(θ), σ', invariantes, J = σ'E). |
| `cap422/plots.py` | Gráficos (setas E/J, elipse, polar, desvio entre E e J). |
| `tests/` | Testes que verificam invariantes, valores numéricos e o caso isotrópico. |
| `outputs/` | Figuras geradas. |
| `docs/CAP422_guia_de_estudo_aula_14-10.pdf` | Guia de estudo: categorias, dual, produto tensorial, tensores (p,q), exercícios com gabarito. |

Rodar os testes: `pytest` (ou o ícone de testes do VS Code).

## Convenção de rotação
O slide usa σ' = R(θ) σ R(θ)ᵀ. Se girar os eixos em vez do tensor, as componentes são Rᵀ σ R (mesmo resultado com θ → −θ). Está implementado em `rotated_sigma(..., convention="slide" | "axes")`.

## Hands-On do Tópico 2: Cálculo Variacional (Prof. Reinaldo)

As três atividades estão na seção "Hands-on Activities" (p. 32-34) da apostila `variational_principle.pdf`, enviada por email em 16/11.

| Arquivo (`topic02_variational/`) | Atividade |
|---|---|
| `handson_1_euler_lagrange.py` / `.ipynb` | Euler-Lagrange simbólico com SymPy (itens 1.1 a 1.6) |
| `handson_2_potential_wells.py` / `.ipynb` | Poços simples e duplo, atrito de Rayleigh, `solve_ivp` (itens 2.1 a 2.7) |
| `handson_3_pendulum.py` / `.ipynb` | Pêndulo, pequenos ângulos, Noether (itens 3.1 a 3.8) |
| `outputs/` | Figuras geradas pelos scripts |

Código reaproveitável novo em `cap422/`: `variational.py` (Euler-Lagrange, Rayleigh) e `dynamics.py` (integração, energia, período do pêndulo).
Os `.ipynb` são gerados a partir dos `.py` com `python tools/py_to_ipynb.py topic02_variational/*.py`.

## Mapa do curso (CAP-422, 2026)

| Tópico | Tema | Professor | Hands-On |
|---|---|---|---|
| 1 | Estruturas matemáticas e tensores | Reinaldo e Heyder | `hands_on_01` (feito) |
| 2 | Cálculo variacional e princípio da mínima ação | Reinaldo | `topic02_variational/` (feito) |
| 3 | Formação de padrões / EDPs não lineares (ordem das aulas pelos emails) | Reinaldo | a fazer |
| 4 | Física de plasmas (aula de apresentação em 27/11) | Reinaldo | a fazer |
| – | Análise complexa e funcional, EDOs e sistemas dinâmicos, perturbação (KAM), grupos de Lie | Heyder | conferir com os slides |

Projeto individual (50% da nota): relatório de 5 a 8 páginas, material computacional e apresentação em 14/12.

## Convenção para os próximos Hands-On
Uma pasta por tópico (`topic03_...`), um script `# %%` por atividade, lógica reaproveitável em `cap422/` e um teste em `tests/` para cada função nova.
