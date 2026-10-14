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

## Próximos Hands-On
Crie uma pasta/arquivo por Hands-On (`hands_on_02.ipynb`, ...) e reaproveite/estenda o pacote `cap422/`.
