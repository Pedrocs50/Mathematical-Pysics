"""Gráficos do Hands-On 1: como J depende da direção de E."""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from .tensors2d import angle_between_deg, current_density, sigma_numeric

E_EXEMPLOS = {"E1=(1,0)": (1.0, 0.0), "E2=(0,1)": (0.0, 1.0), "E3=(1,1)": (1.0, 1.0)}
CORES = ["tab:blue", "tab:orange", "tab:green"]


def _unit(v):
    v = np.asarray(v, dtype=float)
    return v / np.linalg.norm(v)


def plot_vectors_and_ellipse(sx, sy, angles_deg=(0, 30, 45), savepath=None):
    """Linha de cima: setas E (tracejadas) e J (cheias) + elipse da imagem do círculo unitário.
    Linha de baixo: gráfico polar de |J| em função da direção de E (|E| = 1)."""
    n = len(angles_deg)
    fig = plt.figure(figsize=(4.4 * n, 8.6))
    phi = np.linspace(0, 2 * np.pi, 361)
    circle = np.vstack([np.cos(phi), np.sin(phi)])
    lim = 1.25 * max(sx, sy)

    for k, th in enumerate(angles_deg):
        S = sigma_numeric(sx, sy, th)
        ax = fig.add_subplot(2, n, k + 1)
        ellipse = S @ circle
        ax.plot(circle[0], circle[1], color="0.8", lw=1)
        ax.plot(ellipse[0], ellipse[1], color="0.3", lw=1.5, label="σ' · (círculo unitário)")
        for (nome, E), cor in zip(E_EXEMPLOS.items(), CORES):
            E = _unit(E)
            J = current_density(S, E)
            ax.annotate("", xy=E, xytext=(0, 0), arrowprops=dict(arrowstyle="->", color=cor, ls="--", lw=1.2))
            ax.annotate("", xy=J, xytext=(0, 0), arrowprops=dict(arrowstyle="-|>", color=cor, lw=2))
            ax.plot([], [], color=cor, label=f"{nome}  (E tracejado, J cheio)")
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_aspect("equal")
        ax.axhline(0, color="0.9", lw=0.8)
        ax.axvline(0, color="0.9", lw=0.8)
        ax.set_title(f"θ = {th}°")
        if k == 0:
            ax.legend(fontsize=7, loc="lower left")

        axp = fig.add_subplot(2, n, n + k + 1, projection="polar")
        mod = np.linalg.norm(S @ circle, axis=0)
        axp.plot(phi, mod, color="tab:red")
        axp.set_title(f"|J| × direção de E (θ = {th}°)", fontsize=10)
        axp.set_rlim(0, 1.15 * max(sx, sy))

    fig.suptitle(f"Condutividade 2D anisotrópica: σx = {sx}, σy = {sy}", fontsize=13)
    fig.tight_layout()
    if savepath:
        fig.savefig(savepath, dpi=150)
    return fig


def plot_deviation(sx, sy, angles_deg=(0, 30, 45), savepath=None):
    """Ângulo entre J e E em função da direção de E (zero => J paralelo a E)."""
    fig, ax = plt.subplots(figsize=(7, 4))
    phis = np.linspace(0, 360, 721)
    for th in angles_deg:
        S = sigma_numeric(sx, sy, th)
        dev = [angle_between_deg(_unit([np.cos(np.radians(p)), np.sin(np.radians(p))]),
                                 current_density(S, [np.cos(np.radians(p)), np.sin(np.radians(p))]))
               for p in phis]
        ax.plot(phis, dev, label=f"θ = {th}°")
    S_iso = sigma_numeric(sx, sx, 0)
    dev_iso = [angle_between_deg([np.cos(np.radians(p)), np.sin(np.radians(p))],
                                 current_density(S_iso, [np.cos(np.radians(p)), np.sin(np.radians(p))]))
               for p in phis]
    ax.plot(phis, dev_iso, "k--", label=f"isotrópico (σx=σy={sx})")
    ax.set_xlabel("direção de E (graus)")
    ax.set_ylabel("ângulo entre E e J (graus)")
    ax.set_title("Quando J deixa de apontar na direção de E?")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    if savepath:
        fig.savefig(savepath, dpi=150)
    return fig
