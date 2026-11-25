# CAP-422: Física Matemática Computacional (INPE)

Mathematical Physics: A Modern and Computational Approach, com os professores Reinaldo Rosa e Heyder Hey. Hands-On com Python e SymPy, rodando 100% local no VS Code (sem Colab).

## Anotações e teoria

O caderno da disciplina fica em [`docs/`](docs/README.md): linha do tempo, notas de cada aula, teoria e explicação dos Hands-On.

## Atividades

| Tópico (ementa) | Tema | Professor | Pasta | Situação |
|---|---|---|---|---|
| 1 | Estruturas matemáticas e tensores | Reinaldo e Heyder | [`topic01_tensors/`](topic01_tensors/) | implementado |
| 2 | Cálculo variacional e princípio da mínima ação | Reinaldo | [`topic02_variational/`](topic02_variational/) | implementado e [revisado](docs/atividades/topico02_revisao.md); falta reescrever 1.6, 2.7 e 3.7 com as próprias palavras antes de entregar |
| 3 a 6 | Análise complexa e funcional, sistemas dinâmicos, perturbação (KAM), grupos de Lie | Heyder | – | aulas ainda não houve |
| 7 | EDPs não lineares e formação de padrões ("Topic 3" nos emails) | Reinaldo | – | a fazer (falta o enunciado) |
| 8 | Física matemática de plasmas (abertura em 27/11) | Reinaldo | – | a fazer |

A numeração segue a ementa. A situação indica o que está no repositório, não o que já foi entregue.

## Avaliação

| Parte | Peso | Observação |
|---|---|---|
| Presença | 10% | |
| Hands-On | 40% | 8 atividades individuais, 5% cada, uma por tópico |
| Projeto individual | 50% | relatório de 5 a 8 páginas, material computacional e apresentação em **14/12/2026**; tema ainda não escolhido |

## Estrutura

```
docs/                  notas de aula, teoria e explicação das atividades
topic01_tensors/       Tópico 1: guia de estudo, script, notebook e figuras
topic02_variational/   Tópico 2: três scripts, notebooks e figuras
cap422/                código reaproveitável entre os tópicos
tests/                 testes do pacote cap422
tools/py_to_ipynb.py   converte scripts "# %%" em notebooks
```

| Módulo de `cap422/` | Conteúdo | Usado em |
|---|---|---|
| `tensors2d.py` | σ, R(θ), σ', invariantes, J = σ'E | Tópico 1 |
| `plots.py` | gráficos de E e J, elipse, polar, desvio | Tópico 1 |
| `variational.py` | ação, Euler-Lagrange, atrito de Rayleigh | Tópico 2 |
| `dynamics.py` | integração, energia, período do pêndulo | Tópico 2 |

## Como começar

No PowerShell, dentro da pasta do projeto:

```powershell
python -m venv .venv
```

```powershell
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
```

```powershell
.\.venv\Scripts\python.exe -m pytest
```

No Linux ou macOS, o interpretador fica em `.venv/bin/python`.

No VS Code: abra a pasta, instale as extensões recomendadas (Python e Jupyter) e escolha o interpretador `.venv` (Ctrl+Shift+P → *Python: Select Interpreter*).

## Convenção para os próximos Hands-On

- Uma pasta por tópico (`topic03_...`), com um `README.md` que lista o enunciado e a situação.
- Um script com células `# %%` por atividade, com o notebook gerado por `tools/py_to_ipynb.py`. As figuras vão em `outputs/` dentro da pasta do tópico.
- Lógica reaproveitável em `cap422/`, com teste em `tests/` para cada função nova.
- Código e comentários em português.
- Enunciados só da apostila ou dos slides.
