# Tópico 1: estruturas matemáticas e tensores

## Hands-On: tensor de condutividade 2D

| Arquivo | O que é |
|---|---|
| `hands_on_01.ipynb` | Notebook completo, com as respostas das 5 perguntas de discussão. |
| `hands_on_01.py` | Os cálculos e as figuras em script com células `# %%`. Não contém as respostas das perguntas. |
| `outputs/` | Figuras geradas. |
| `CAP422_guia_de_estudo_aula_14-10.pdf` | Guia de estudo: categorias, dual, produto tensorial, tensores (p,q), exercícios com gabarito. |

Explicação dos resultados: [como explicar](../docs/atividades/topico01_como_explicar.md).

Código usado: `cap422/tensors2d.py` e `cap422/plots.py`. Testes: `tests/test_tensors2d.py`.

Para rodar tudo, a partir da raiz do repositório:

```powershell
.\.venv\Scripts\python.exe topic01_tensors/hands_on_01.py
```

> **Atenção:** neste tópico, não gere o notebook com `tools/py_to_ipynb.py`. Ele sobrescreveria o `hands_on_01.ipynb` e apagaria as respostas e as saídas, que só existem no notebook.

## Convenção de rotação

O slide usa σ' = R(θ) σ R(θ)ᵀ. Se girar os eixos em vez do tensor, as componentes são Rᵀ σ R (mesmo resultado com θ → −θ). Está implementado em `rotated_sigma(..., convention="slide" | "axes")`.
