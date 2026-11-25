# Teoria 3: formação de padrões

[Índice](../README.md) · [Aula de 25/11](../aulas/2026-11-25_padroes_parte_A.md) · [Anterior](02_calculo_variacional.md)

> **Resumo.** Sistemas extensos fora do equilíbrio criam estruturas (células, espirais, turbulência). Estudamos como **representar** esses dados, como **gerá-los** com modelos e como **medi-los**.

Fontes: os três artigos enviados para o "Topic 3" (na ementa, tópico 7). Trechos marcados como **Complemento** são explicação de apoio.

**Sumário:** [Panorama](#1-panorama) · [GNL](#2-redes-numéricas-generalizadas-31) · [Mapas acoplados](#3-redes-de-mapas-acoplados-31) · [GPA](#4-descritores-espaço-temporais-32) · [Ginzburg–Landau](#5-modelo-de-edp-ginzburglandau-33) · [Navier–Stokes](#6-modelo-de-edp-navierstokes-33)

---

## 1. Panorama

> **Complemento.**

Três formas de gerar um padrão no computador:

| Modelo | Espaço | Tempo | Variável | Exemplo aqui |
|---|---|---|---|---|
| rede de mapas acoplados | discreto | discreto | contínua | seção 3 |
| EDP de amplitude | contínuo | contínuo | complexa | Ginzburg–Landau |
| EDP de fluido | contínuo | contínuo | vetorial | Navier–Stokes |

Em todos, o padrão nasce da disputa entre uma **dinâmica local** (não linear) e um **transporte** que conecta vizinhos.

---

## 2. Redes numéricas generalizadas (3.1)

*Veronese et al., 2009.*

**Problema:** um mesmo fenômeno (uma explosão solar, por exemplo) é medido por muitos instrumentos, com dimensões e resoluções diferentes. Como organizar tudo?

**Proposta:** descrever cada série por três tipos de coeficientes.

$$\pounds=f(\kappa,\lambda,\mu_p)$$

| Símbolo | Nome | Responde |
|---|---|---|
| κ | grau variacional | de quantas variáveis os dados dependem? |
| λ | extensões | quantos pontos há em cada uma? |
| μp | parâmetros pós-analíticos | o que a análise mediu nesses dados? |

### O grau variacional κ

| κ | Variáveis | Exemplo |
|---|---|---|
| 1 | t | só o tempo |
| 2 | U(t) | série temporal |
| 3 | U(t, x) | perfil que evolui |
| 4 | U(t, x, y) | sequência de imagens |
| 5 | U(t, x, y, z) | sequência de volumes |
| ≥ 6 | variáveis acopladas extras | genômica, SIG, cosmologia |

### Como ler a notação

| Escrita | Significado |
|---|---|
| £₂, λ₀ = 10⁴ | série temporal com 10 000 pontos |
| £₂,₁₀⁴ (0,28; −1,66) | a mesma, com autocorrelação μ₁ = 0,28 e índice espectral μ₂ = −1,66 |
| £₄,₁₀²,₆₄,₆₄ | 100 imagens de 64 × 64 |
| £₅,₁₀²,₆₄,₆₄,₆₄ | 100 volumes de 64 × 64 × 64 |

### Os parâmetros μ

- **μ₁, autocorrelação** (entre −1 e 1): detecta repetição; perto de 0 para ruído branco.
- **μ₂, índice do espectro de potência**: mede o decaimento do tipo lei de potência (ruído 1/f).
- Para imagens (κ = 4): funções de correlação espacial, funcionais de Minkowski e momentos de gradiente (GPA, seção 4).
- Outros citados: expoente de Hurst, entropia de Kolmogorov–Sinai, dimensões fractais.

### Data System Representation

Uma grade com todas as GNLs de um sistema: **linhas por λ₀, colunas por κ**, com totais nas margens. O estudo de caso é a atividade solar de 6 de junho de 2000, com dados do Sol, do meio interplanetário e da magnetosfera. O artigo apresenta o algoritmo como trabalho em desenvolvimento.

---

## 3. Redes de mapas acoplados (3.1)

> **Complemento.** O tema está na agenda da aula, mas os artigos enviados só o citam de passagem. Confira a definição usada pelo professor.

Cada sítio da rede tem um valor que evolui por um mapa não linear f e troca informação com os vizinhos. A forma difusiva clássica (Kaneko), em uma dimensão:

$$x_{t+1}(i)=(1-\varepsilon)\,f\big(x_t(i)\big)+\frac{\varepsilon}{2}\Big[f\big(x_t(i-1)\big)+f\big(x_t(i+1)\big)\Big]$$

| Peça | Papel |
|---|---|
| f, por exemplo f(x) = r x (1 − x) | dinâmica local, que pode ser caótica |
| ε entre 0 e 1 | força do acoplamento |
| bordas periódicas | a rede "dá a volta" |

Em duas dimensões, a média usa os 4 vizinhos, com peso ε/4. Variando r e ε surgem regimes ordenados, intermitentes e turbulentos. O resultado é um dado do tipo £₃ (rede 1D) ou £₄ (rede 2D).

---

## 4. Descritores espaço-temporais (3.2)

**Análise de Padrões de Gradiente (GPA).** Em vez de olhar o histograma dos valores, a GPA olha **onde o campo de gradiente quebra a simetria**.

### O algoritmo

```
matriz M  →  gradiente ∇M  →  remove pares simétricos  →  rede assimétrica  →  momentos G1 a G4
```

1. Calcula o gradiente por diferenças centrais (de primeira ordem nas bordas).
2. Remove pares de vetores **simétricos**: mesmo módulo, fases opostas, mesma distância do centro.
3. Sobra a **rede assimétrica**, com V_A dos V vetores.
4. Calcula os momentos.

### Os momentos

| Momento | Mede | Fórmula |
|---|---|---|
| G1 | geometria (triangulação de Delaunay, T_A conexões) | (T_A − V_A) / V_A |
| G2 | módulos e alinhamento | ver abaixo |
| G3 | fases | ver abaixo |
| G4 | entropia complexa | −Σ zᵢ ln zᵢ, com zᵢ = \|vᵢ\| e^(iφᵢ) |

$$G_2=\frac{V_A}{V}\left(1-\frac{\lvert\sum_i\mathbf v_i\rvert}{2\sum_i\lvert\mathbf v_i\rvert}\right)$$

| Padrão | G2 |
|---|---|
| todo simétrico | 0 |
| laminar (todos assimétricos e alinhados) | 0,5 |
| assimétricos apontando para todo lado | V_A / V (perto de 1 se quase todos são assimétricos) |

**G3 refinado** (a novidade do artigo de 2023, equação 19):

$$G_3=\frac12\left(\frac{V_A}{V}+\frac{1}{2V_A}\sum_i\big(\mathbf u_i\cdot\mathbf v_i+1\big)\right)$$

Aqui uᵢ é o vetor unitário na direção do ponto vista do centro de simetria, e vᵢ é o vetor unitário na direção da fase. O produto vale 1 se estão alinhados, 0 se são perpendiculares e −1 se são opostos. Isso evita o salto entre os ângulos 0 e 2π da versão antiga.

Segundo o artigo, G2 e G3 distinguem melhor flutuações finas do que G1 e G4 e se complementam.

---

## 5. Modelo de EDP: Ginzburg–Landau (3.3)

*Sautter, Rosa e Pontes, 2023.* Prévia da parte B.

A equação complexa de Ginzburg–Landau estocástica, em duas versões:

$$\partial_tA=(1+ib)\nabla^2A+A-(1+ic)\lvert A\rvert^2A+\sigma\,\partial_t\eta_\beta\qquad\text{(ruído aditivo)}$$

$$\partial_tA=(1+ib)\nabla^2A+A-(1+ic)\lvert A\rvert^2A+\sigma A\,\partial_t\eta_\beta\qquad\text{(ruído multiplicativo)}$$

| Termo | Papel |
|---|---|
| (1 + ib)∇²A | difusão e dispersão no espaço |
| A | crescimento linear |
| −(1 + ic)\|A\|²A | saturação não linear |
| ruído η_β | flutuação com espectro 1/f^β |

| Ruído | β |
|---|---|
| branco | 0 |
| rosa | 1 |
| vermelho (browniano) | 2 |

**Numérico:** método pseudo-espectral de Fourier no espaço, Runge–Kutta–Fehlberg adaptativo no tempo, bordas periódicas. Caso estudado: b = 0, c = 1,5, σ = 100, β = 2, Δt = 0,05, Δh = 1, tolerância 10⁻⁴, grade 128 × 128.

### Três regimes (sem ruído)

| Regime | Instante mostrado | O que se vê |
|---|---|---|
| acreção | T = 2 | formam-se células |
| formação de espirais | T = 12 | defeitos espirais |
| relaxação | T = 57 | as espirais se acomodam |

### Resultados

- Uma análise de Lyapunov da parte reativa indica que **só o ruído vermelho** produz instabilidade relevante.
- O ruído **multiplicativo** preserva as grandes estruturas espirais, parecido com o caso sem ruído. O **aditivo** inibe a formação delas.
- No plano G2 × G3, os três modelos (sem ruído, aditivo, multiplicativo) ficam em regiões separáveis.

---

## 6. Modelo de EDP: Navier–Stokes (3.3)

*Mendes e Rosa, ENMC 2025.* Prévia da parte B.

**DNS** (simulação numérica direta) resolve todas as escalas do escoamento, sem modelo de turbulência. É a abordagem mais fiel e a mais cara (comparada com RANS e LES).

$$\frac{\partial\mathbf u}{\partial t}-\mathbf u\times\boldsymbol\omega=-\nabla P+\nu\nabla^2\mathbf u+\mathbf f,\qquad\nabla\cdot\mathbf u=0$$

com ω = ∇ × u e P = p + ½ u·u, em um cubo periódico [0, 2π]³.

### A receita pseudo-espectral

| Passo | Por quê |
|---|---|
| derivar no espaço de Fourier | derivada vira multiplicação por ik, com alta precisão |
| calcular o produto não linear no espaço físico | evita uma convolução cara |
| regra dos 2/3 | remove o *aliasing* |
| projeção solenoidal | garante ∇·u = 0 |
| fator integrante | trata a viscosidade de forma exata |
| Runge–Kutta de 4ª ordem | avança no tempo |

Implementação em Python com mpi4py, CuPy e NCCL (várias GPUs).

### Validação: vórtice de Taylor–Green, Re = 1600, contra dados da NASA (512³)

| Grade | Erro na energia cinética | Erro na enstrofia | Erro na dissipação |
|---|---|---|---|
| 64³ | 12,77% | 10,88% | 10,89% |
| 128³ | 3,21% | 4,76% | 4,76% |
| 256³ | 0,83% | 2,53% | 2,53% |
| 512³ | 0,02% | 0,11% | 0,11% |

### Desempenho

- Com 8 GPUs V100, o tempo caiu de 17 h 33 min 24 s (um núcleo de CPU) para 6 min 03 s, cerca de **174 vezes**.
- Uma GPU P100 equivaleu a cerca de 49 núcleos; uma V100-SXM2, a cerca de 84.
- Em turbulência forçada (256³, Re = 1800, forçamento em k = 4 a 6), o espectro seguiu **k^(−5/3)**.

---

## Teste-se

<details>
<summary>Qual o κ de um filme de magnetogramas do Sol? Como escrever 100 quadros de 256 × 256?</summary>

κ = 4: £₄ com λ₀ = 100, λ₁ = λ₂ = 256.
</details>

<details>
<summary>Por que a GPA remove os pares simétricos?</summary>

Porque o que ela quer medir é a quebra de simetria. Um padrão totalmente simétrico fica sem vetores e tem G2 = 0.
</details>

<details>
<summary>Qual a diferença entre ruído aditivo e multiplicativo?</summary>

O aditivo soma a perturbação independentemente do estado. No multiplicativo ela é proporcional a A, então some onde A é zero.
</details>
