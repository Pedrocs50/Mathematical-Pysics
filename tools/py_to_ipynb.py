"""Converte um script com células '# %%' (e '# %% [markdown]') em notebook .ipynb.

Uso:  python tools/py_to_ipynb.py topic02_variational/handson_1_euler_lagrange.py
Gera o .ipynb ao lado do .py (sem saídas; abra no VS Code e rode as células).
"""
from __future__ import annotations

import sys
from pathlib import Path

import nbformat as nbf


def convert(path: Path) -> Path:
    cells, kind, buf = [], None, []

    def flush():
        if kind is None or not "".join(buf).strip():
            return
        text = "".join(buf).strip("\n")
        if kind == "markdown":
            lines = [ln[2:] if ln.startswith("# ") else ln.lstrip("#") for ln in text.splitlines()]
            cells.append(nbf.v4.new_markdown_cell("\n".join(lines)))
        else:
            cells.append(nbf.v4.new_code_cell(text))

    for line in path.read_text(encoding="utf-8").splitlines(keepends=True):
        if line.startswith("# %%"):
            flush()
            kind, buf = ("markdown" if "[markdown]" in line else "code"), []
        else:
            buf.append(line)
    flush()
    nb = nbf.v4.new_notebook(cells=cells)
    nb.metadata["kernelspec"] = {"display_name": "Python 3", "language": "python", "name": "python3"}
    out = path.with_suffix(".ipynb")
    nbf.write(nb, out)
    return out


if __name__ == "__main__":
    for p in sys.argv[1:]:
        print("gerado:", convert(Path(p)))
