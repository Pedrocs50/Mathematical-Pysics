# %% [markdown]
# # Hands-On 1: Tensor de condutividade 2D (CAP-422)
# Rode célula por célula no VS Code (botão "Run Cell" acima de cada `# %%`)
# ou tudo de uma vez: `python hands_on_01.py`.

# %%
import matplotlib

if __name__ == "__main__":
    matplotlib.use("Agg")  # sem abrir janela quando roda como script
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp

from cap422 import tensors2d as t
from cap422.plots import plot_deviation, plot_vectors_and_ellipse

sp.init_printing(use_unicode=True)

# %% 1) Definir sigma simbolicamente
sigma = t.sigma_symbolic()
print("1) sigma =")
sp.pprint(sigma)

# %% 2) Matriz de rotação e sigma'(theta) = R sigma R^T
R = t.rotation_matrix()
sigma_p = t.rotated_sigma(sigma)
print("\n2) R(theta) =")
sp.pprint(R)
print("\n   sigma'(theta) = R sigma R^T =")
sp.pprint(sigma_p)

# %% 3) Simplificação: forma de ângulo duplo
sigma_p2 = t.rotated_sigma_double_angle()
print("\n3) Forma de ângulo duplo:")
sp.pprint(sigma_p2)
print("   Iguais às anteriores?", sp.simplify(sigma_p - sigma_p2) == sp.zeros(2))

# %% 4) Invariantes
d_tr, d_det = t.invariants_difference(sigma, sigma_p)
print("\n4) tr(sigma') - tr(sigma) =", d_tr, "| det(sigma') - det(sigma) =", d_det)
print("   tr  =", sp.simplify(sigma_p.trace()), "| det =", sp.simplify(sigma_p.det()))

# %% 5) Valores numéricos (sigma_x = 2, sigma_y = 1)
SX, SY = 2.0, 1.0
print("\n5) sigma' numérico para sigma_x=2, sigma_y=1:")
for th in (0, 30, 45):
    S = t.sigma_numeric(SX, SY, th)
    print(f"   theta = {th:>2}°  ->\n{np.round(S, 4)}  tr={np.trace(S):.4f}  det={np.linalg.det(S):.4f}")

# %% 6) J = sigma' E com E = (1, 1)
E = np.array([1.0, 1.0])
print("\n6) E = (1, 1)")
for th in (0, 30, 45):
    S = t.sigma_numeric(SX, SY, th)
    J = t.current_density(S, E)
    ang = t.angle_between_deg(E, J)
    print(f"   theta = {th:>2}°  J = {np.round(J, 4)}  ângulo(E -> J) = {ang:+.2f}°")

# %% 7) Visualização
fig1 = plot_vectors_and_ellipse(SX, SY, (0, 30, 45), savepath="outputs/J_vs_E.png")
fig2 = plot_deviation(SX, SY, (0, 30, 45), savepath="outputs/desvio_E_J.png")
print("\n7) Figuras salvas em outputs/")
if __name__ != "__main__":
    plt.show()
