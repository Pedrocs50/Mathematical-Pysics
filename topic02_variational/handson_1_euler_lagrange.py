# %% [markdown]
# # Tópico 2, Hands-On 1: Equações de Euler-Lagrange simbólicas (CAP-422)
# Objetivo: obter a equação de movimento a partir de uma Lagrangiana, usando derivação simbólica.
# Rode célula por célula no VS Code ou tudo de uma vez: `python topic02_variational/handson_1_euler_lagrange.py`

# %%
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1] if "__file__" in globals() else Path.cwd().parent
sys.path.insert(0, str(ROOT))

import sympy as sp
from sympy.calculus.euler import euler_equations

from cap422.variational import action, euler_lagrange

sp.init_printing(use_unicode=True)

# %% [markdown]
# ## 1.1 Símbolos
# `x` é a variável independente, `y(x)` é uma função e `y'(x)` a derivada dela.

# %%
x = sp.Symbol("x", real=True)
y = sp.Function("y")(x)
yp = sp.diff(y, x)
print("1.1) y =", y, "| y' =", yp)

# %% [markdown]
# ## 1.2 A ação
# S[y] = ∫ L(x, y, y') dx.
# **y(x) é uma função**: recebe um número x e devolve um número.
# **S[y] é um funcional**: recebe a função inteira y e devolve um único número.
# Por isso o cálculo variacional procura a *função* y que deixa S estacionária.

# %%
L = sp.Function("L")(x, y, yp)  # L genérica, só para mostrar a forma
x0, x1 = sp.symbols("x0 x1")
print("1.2) S[y] =")
sp.pprint(action(L, x, x0, x1))

# %% [markdown]
# ## 1.3 Rotina reutilizável
# A rotina está em `cap422/variational.py::euler_lagrange(L, y, x)` e calcula d/dx(∂L/∂y') − ∂L/∂y.
# Conferimos contra a função pronta do SymPy (`euler_equations`) em cada teste abaixo.

# %%
def resolver(L, y, x):
    """Aplica a nossa rotina, confere com o SymPy e devolve a equação diferencial."""
    meu = euler_lagrange(L, y, x)
    ref = euler_equations(L, y, x)[0].lhs  # o SymPy usa o sinal oposto: ∂L/∂y − d/dx ∂L/∂y'
    assert sp.simplify(meu + ref) == 0, "diferente do euler_equations do SymPy"
    return sp.Eq(meu, 0)

# %% [markdown]
# ## 1.4 L = ½ y'² + ½ ω² y²   (atenção ao sinal de +)

# %%
w = sp.Symbol("omega", positive=True)
L1 = sp.Rational(1, 2) * yp**2 + sp.Rational(1, 2) * w**2 * y**2
eq1 = resolver(L1, y, x)
print("1.4) Equação de Euler-Lagrange:")
sp.pprint(eq1)
sol1 = sp.dsolve(eq1, y)
print("     Solução simbólica:")
sp.pprint(sol1)

# %% [markdown]
# ## 1.5 L = ½ m y'² − ½ k y²

# %%
m, k = sp.symbols("m k", positive=True)
L2 = sp.Rational(1, 2) * m * yp**2 - sp.Rational(1, 2) * k * y**2
eq2 = resolver(L2, y, x)
print("1.5) Equação de Euler-Lagrange:")
sp.pprint(eq2)
sol2 = sp.dsolve(eq2, y)
print("     Solução simbólica:")
sp.pprint(sol2)

# %% [markdown]
# ### Comparação dos dois sistemas
# - Em 1.4 a equação é **y'' − ω² y = 0**. As soluções são exponenciais reais (cosh e senh):
#   quanto mais longe de y = 0, mais a força empurra para longe. O ponto y = 0 é **instável**.
#   Isso acontece porque o sinal de + na Lagrangiana dá um "potencial" −½ω²y², um morro invertido.
# - Em 1.5 a equação é **m y'' + k y = 0**. As soluções são senos e cossenos com frequência
#   √(k/m): o sistema oscila em torno de y = 0, que é **estável** (oscilador harmônico).
# - Os dois têm a mesma forma, y'' = (constante) · y. A diferença está só no sinal da constante.

# %% [markdown]
# ## 1.6 Euler-Lagrange e a primeira variação
# Compare a curva y com uma curva vizinha y + εη, onde η(x0) = η(x1) = 0 (as pontas ficam fixas).
# Para cada η escolhido, S[y + εη] é uma função comum do número ε. A **primeira variação** é a
# derivada dessa função em ε = 0:
#
# $$\delta S=\frac{d}{d\varepsilon}S[y+\varepsilon\eta]\Big|_{\varepsilon=0}
# =\int_{x_0}^{x_1}\left(\frac{\partial L}{\partial y}\,\eta+\frac{\partial L}{\partial y'}\,\eta'\right)dx$$
#
# Integrando por partes o segundo termo:
#
# $$\delta S=\left[\frac{\partial L}{\partial y'}\,\eta\right]_{x_0}^{x_1}
# +\int_{x_0}^{x_1}\left(\frac{\partial L}{\partial y}-\frac{d}{dx}\frac{\partial L}{\partial y'}\right)\eta\,dx$$
#
# O primeiro termo (o de borda) é zero porque η se anula nas duas pontas. Para y ser estacionária,
# δS = 0 tem de valer para **qualquer** η. Pelo lema fundamental do cálculo variacional, isso só é
# possível se o parêntese for zero em todo ponto do intervalo. Com o sinal trocado, esse parêntese
# é a expressão do item 1.3, então d/dx(∂L/∂y') − ∂L/∂y = 0.
# Ou seja: **a equação de Euler-Lagrange é a condição de primeira variação nula.** Ela garante que
# S é estacionária, não que seja mínima.
