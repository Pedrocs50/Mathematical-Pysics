# Teoria 1: estruturas matemáticas e tensores

[Índice](../README.md) · [Aula de 14/10](../aulas/2026-10-14_tensores.md) · [Hands-On](../atividades/topico01_como_explicar.md) · [Próxima teoria](02_calculo_variacional.md)

> **Resumo.** Categorias dão a linguagem (objetos e setas). Espaços vetoriais dão o terreno. Tensores são objetos com uma lei de transformação bem definida quando a base muda.

Fonte: [guia de estudo](../../topic01_tensors/CAP422_guia_de_estudo_aula_14-10.pdf), que tem mais exemplos e 8 exercícios com gabarito.

**Sumário:** [Categorias](#1-categorias) · [Escada de espaços](#2-a-escada-de-espaços) · [Base e dual](#3-vetores-base-e-dual) · [Produto tensorial](#4-produto-tensorial) · [Tensores](#5-tensores-do-tipo-pq) · [Condutividade](#6-aplicação-o-tensor-de-condutividade)

---

## 1. Categorias

Uma **categoria** tem quatro ingredientes:

1. objetos;
2. morfismos (setas) entre objetos;
3. composição: dadas f: A → B e g: B → C, existe g ∘ f: A → C;
4. identidade em cada objeto.

E duas leis: a composição é associativa, e a identidade não altera nada.

| Categoria | Objetos | Morfismos | Isomorfismo |
|---|---|---|---|
| Set | conjuntos | funções | bijeção |
| Grp | grupos | homomorfismos | isomorfismo de grupos |
| Vect | espaços vetoriais | aplicações lineares | bijeção linear |
| Top | espaços topológicos | funções contínuas | homeomorfismo |
| Met | espaços métricos | isometrias ou contínuas, conforme a convenção | isometria bijetiva |

- **Invariante:** o que todo isomorfismo preserva. Gênero em Top; dimensão em Vect; traço, determinante e autovalores de um operador.
- **Xícara e rosca:** ambas têm gênero 1, então são "a mesma coisa" em Top, mas não em Met.
- **Funtor:** leva objetos e setas de uma categoria para outra, respeitando composição e identidade. Se inverte as setas, é **contravariante**.
- **Transformação natural:** uma "tradução" entre funtores que não depende de escolhas. A aplicação V → V\*\* é natural; um isomorfismo V → V\* exige escolher uma base ou um produto interno.

## 2. A escada de espaços

| Espaço | Estrutura acrescentada | O que passa a ser possível |
|---|---|---|
| topológico | conjuntos abertos | continuidade |
| métrico | distância | medir separação |
| normado | vetores e norma | comprimento |
| com produto interno | forma simétrica positiva | ângulos |
| Hilbert | completude | limites, bases de Fourier |

A álgebra tensorial só precisa do espaço vetorial. Métrica e produto interno são estrutura **extra**.

## 3. Vetores, base e dual

Escolhida uma base {eᵢ}, todo vetor se escreve v = vⁱ eᵢ (índice repetido = soma).

> As componentes dependem da base. O vetor, não.

Se a nova base é e' = e P, então:

| Objeto | Componentes novas | Nome |
|---|---|---|
| vetor | v' = P⁻¹ v | contravariante (varia ao contrário da base) |
| operador linear | A' = P⁻¹ A P | similaridade |
| covetor (linha ω) | ω' = ω P | covariante (varia com a base) |

O **dual** V\* é o conjunto das funções lineares ω: V → ℝ. A base dual obedece εⁱ(eⱼ) = δⁱⱼ.

O número ω(v) = ωᵢ vⁱ não depende da base: P e P⁻¹ se cancelam. Essa soma de um índice de cima com um de baixo é a **contração**.

**Exemplo físico:** o gradiente df é um covetor; aplicado a v, dá a derivada de f na direção de v.

## 4. Produto tensorial

V ⊗ W é gerado pelos símbolos v ⊗ w, lineares em cada entrada.

- Base: eᵢ ⊗ fⱼ. Dimensão: dim V × dim W.
- Elementos v ⊗ w equivalem a matrizes de **posto 1**. Nem todo tensor é assim.
- **Propriedade universal:** toda aplicação bilinear V × W → U equivale a uma única aplicação **linear** V ⊗ W → U.
- Em dimensão finita, operadores V → V são o mesmo que elementos de V ⊗ V\*.

## 5. Tensores do tipo (p,q)

$$T\in\underbrace{V\otimes\cdots\otimes V}_{p}\otimes\underbrace{V^{*}\otimes\cdots\otimes V^{*}}_{q}$$

A ordem é p + q, e o número de componentes é n^(p+q).

| Tipo | Objeto | Exemplo |
|---|---|---|
| (0,0) | escalar | temperatura, massa, potencial elétrico |
| (1,0) | vetor | velocidade |
| (0,1) | covetor | gradiente, campo elétrico |
| (1,1) | operador | transformação linear, condutividade (classificação do slide) |
| (0,2) | forma bilinear | métrica |
| (2,0) | contravariante de ordem 2 | tensão |

Lei de transformação em forma matricial (nova base e' = e P):

| Tipo | Lei |
|---|---|
| (1,1) | T' = P⁻¹ T P |
| (0,2) | T' = Pᵀ T P |
| (2,0) | T' = P⁻¹ T P⁻ᵀ |

Para rotações, P⁻¹ = Pᵀ, e as três leis coincidem.

> **Detalhe que vale saber.** O slide trata a condutividade como (1,1). Como J é vetor e E é naturalmente covetor, a rigor σ seria (2,0). Em base ortonormal os números e a lei de rotação são os mesmos.

## 6. Aplicação: o tensor de condutividade

$$\mathbf J=\sigma\,\mathbf E,\qquad\sigma'=R(\theta)\,\sigma\,R(\theta)^{T}$$

Com σ = diag(σx, σy), o ângulo duplo separa as duas partes:

$$\sigma'=\frac{\sigma_x+\sigma_y}{2}\,I+\frac{\sigma_x-\sigma_y}{2}
\begin{pmatrix}\cos2\theta&\sin2\theta\\ \sin2\theta&-\cos2\theta\end{pmatrix}$$

A primeira parte é isotrópica e não gira. A segunda mede a anisotropia e gira com 2θ. Veja o desenvolvimento em [como explicar o Hands-On](../atividades/topico01_como_explicar.md).

## Teste-se

<details>
<summary>O que faz de um array um tensor?</summary>

A lei de transformação. O objeto é uma aplicação multilinear; o array é a representação dele numa base.
</details>

<details>
<summary>Dim V = 2 e dim W = 3. Qual a dimensão de V ⊗ W? E de V ⊗ V*?</summary>

6 e 4.
</details>

<details>
<summary>Dê um funtor covariante e um contravariante de Vect em Vect.</summary>

Covariante: V ⊗ (−), ou a identidade. Contravariante: o dual, V ↦ V*.
</details>
