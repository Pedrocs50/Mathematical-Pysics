# Teoria 2: cálculo variacional e mínima ação

[Índice](../README.md) · [Aula de 16/11](../aulas/2026-11-16_calculo_variacional.md) · [Hands-On](../atividades/topico02_como_explicar.md) · [Anterior](01_tensores.md) · [Próxima](03_formacao_de_padroes.md)

> **Resumo.** Um funcional dá um número para cada função. Exigir que sua primeira variação seja zero leva à equação de Euler–Lagrange. Quando o funcional é a ação S = ∫(T − V) dt, essa equação é a 2ª lei de Newton.

Fonte: apostila *The Variational Principle* (R. R. Rosa, novembro de 2026) e slides do Tópico 2.

**Sumário:** [Funções](#1-começando-por-funções) · [Funcionais](#2-funcionais) · [Euler–Lagrange](#3-a-equação-de-eulerlagrange) · [Corda](#4-o-exemplo-da-corda) · [Primeira integral](#5-primeira-integral) · [Mecânica](#6-mecânica-o-princípio-da-ação-estacionária) · [Noether](#7-noether-simetria-e-conservação) · [Mínimo?](#8-estacionário-ou-mínimo) · [Atrito](#9-atrito-de-rayleigh)

---

## 1. Começando por funções

No mínimo ou máximo x_m de f, a variação é nula em primeira ordem:

$$\delta f=f(x_m+\delta x)-f(x_m)=f'(x_m)\,\delta x+O(\delta x^2)$$

Pedir δf = 0 para **qualquer** δx é o mesmo que f'(x_m) = 0.

Com várias variáveis, a expansão de Taylor fica compacta:

$$f(\mathbf r+\mathbf a)=e^{\,\mathbf a\cdot\nabla}f(\mathbf r)=f+\mathbf a\cdot\nabla f+\tfrac12\,\mathbf a^{T}H\,\mathbf a+\dots$$

onde H é a **Hessiana** (a matriz das segundas derivadas).

| f'(x) = 0 e | Tipo de ponto |
|---|---|
| f'' > 0 | mínimo |
| f'' < 0 | máximo |
| f'' = 0 | o teste não decide |

## 2. Funcionais

Um **funcional** recebe uma função inteira e devolve um número. Pense nele como uma função de infinitas variáveis, uma para cada ponto x.

| Notação | O que é |
|---|---|
| f(x) | função: número entra, número sai |
| I[f] | funcional: função entra, número sai |

Exemplo da apostila: a energia potencial de uma corda de densidade σ com formato f(x).

$$U[f]=\int_{x_a}^{x_b}\sigma g\,f(x)\sqrt{1+f'(x)^2}\,dx$$

### δ não é d

```
 y
 |                      ● y + δy   (curva variada)
 |                 .·´
 |            ●··´      ↕ δy : compara duas curvas no mesmo x
 |        __--○─────────  y(x)     (curva original)
 |   __--          dy : anda sobre a mesma curva, de x para x + dx
 └──────────┼──────┼──────── x
            x    x+dx
```

## 3. A equação de Euler–Lagrange

Queremos a função y que torna estacionário

$$S[y]=\int_{x_0}^{x_1}L(x,y,y')\,dx$$

**Passo 1.** Perturbe: y → y + εη, com η(x₀) = η(x₁) = 0.

**Passo 2.** Defina a primeira variação como uma derivada comum em ε:

$$\delta S=\frac{d}{d\varepsilon}S[y+\varepsilon\eta]\Big|_{\varepsilon=0}
=\int_{x_0}^{x_1}\left(\frac{\partial L}{\partial y}\,\eta+\frac{\partial L}{\partial y'}\,\eta'\right)dx$$

**Passo 3.** Integre o segundo termo por partes:

$$\delta S=\left[\frac{\partial L}{\partial y'}\,\eta\right]_{x_0}^{x_1}
+\int_{x_0}^{x_1}\left(\frac{\partial L}{\partial y}-\frac{d}{dx}\frac{\partial L}{\partial y'}\right)\eta\,dx$$

**Passo 4.** O termo de borda é zero porque η se anula nas pontas.

**Passo 5.** δS = 0 para **qualquer** η. Pelo lema fundamental do cálculo variacional, o parêntese é zero em todo ponto:

$$\boxed{\frac{d}{dx}\frac{\partial L}{\partial y'}-\frac{\partial L}{\partial y}=0}$$

Em geral é uma equação diferencial de segunda ordem. Com várias funções, há uma equação para cada uma.

## 4. O exemplo da corda

A apostila usa a corda pendurada para mostrar que **a pergunta precisa estar certa**.

**Pergunta incompleta:** qual formato minimiza a energia? A resposta que sai é uma corda esticada indefinidamente para baixo (|f'| → ∞). A matemática está certa: nada fixou o comprimento.

**Pergunta correta:** qual formato minimiza a energia **com comprimento l₀**?

### Multiplicadores de Lagrange

Para extremizar S com a restrição φ = 0, extremize S − λφ sem restrição.

| Caso | Resultado |
|---|---|
| retângulo de maior área com x + y = a | ∇S = λ∇φ dá x = y = a/2, área a²/4 |
| corda com comprimento fixo | L = σg (f − ζ) √(1 + f'²), com ζ = λ/σg |

A solução é a **catenária**:

$$f(x)=\zeta+\sqrt{C_2}\,\cosh\!\left(\frac{x-x_0}{\sqrt{C_2}}\right)$$

As três constantes vêm só da geometria (as duas pontas e o comprimento). A densidade não aparece: corrente de aço e fio de algodão tomam a mesma forma.

A apostila ainda trata a **restrição local** (uma condição em cada ponto, com λ(l)) e a **mola pendurada**, cuja forma vira a catenária quando a mola fica infinitamente rígida.

## 5. Primeira integral

Se L **não depende explicitamente de x**, então ao longo das soluções:

$$H=y'\,\frac{\partial L}{\partial y'}-L=\text{constante}$$

Prova em uma linha:

$$\frac{dH}{dx}=y'\left(\frac{d}{dx}\frac{\partial L}{\partial y'}-\frac{\partial L}{\partial y}\right)=0$$

Isso baixa a ordem da equação e é o caso mais simples do teorema de Noether.

## 6. Mecânica: o princípio da ação estacionária

Troque x por t, y por q e escreva L = T − V:

$$S[q]=\int_{t_1}^{t_2}L(q,\dot q,t)\,dt,\qquad\delta S=0$$

Para uma partícula, L = ½ m ẋ² − V(x):

$$m\ddot x=-\frac{dV}{dx}=F$$

| Derivada de L | Significado |
|---|---|
| ∂L/∂q̇ | momento conjugado p |
| ∂L/∂q | força generalizada |

Euler–Lagrange fica com a cara de Newton: dp/dt = ∂L/∂q.

### Hamilton e Hamilton–Jacobi

$$H=\sum_i p_i\dot q_i-L,\qquad p_i=\frac{\partial L}{\partial\dot q_i}$$

Para a partícula, H = p²/2m + V = T + V. A formulação de Hamilton–Jacobi resume tudo em uma EDP para a função principal S:

$$\frac{\partial S}{\partial t}+H\!\left(q,\frac{\partial S}{\partial q},t\right)=0$$

| Grandeza | Expressão | Unidade |
|---|---|---|
| T, V, L, H | ½mẋ², V(x), T − V, T + V | J |
| ação S | ∫L dt | J·s |

### Para que serve, se Newton já resolve?

- L é um escalar, então trocar de coordenadas é fácil.
- Restrições entram na escolha das coordenadas.
- Leis de conservação aparecem direto da forma de L.
- A mesma linguagem se estende a meios contínuos, teoria de campos e mecânica quântica.

**Óptica:** o princípio de Fermat (tempo mínimo) é outro problema variacional, T[f] = ∫ (n/c) √(1 + f'²) dx.

## 7. Noether: simetria e conservação

> Toda simetria contínua da ação implica uma lei de conservação.

| Simetria | Transformação | Conserva |
|---|---|---|
| homogeneidade do tempo | t → t + ε | energia |
| homogeneidade do espaço | xⁱ → xⁱ + εⁱ | momento linear |
| isotropia do espaço | xⁱ → Rⁱⱼ xʲ | momento angular |

## 8. Estacionário ou mínimo?

δS = 0 só garante que a ação é **estacionária**. Para saber o tipo, olhe a segunda variação. Com L = ½ m ẋ² − V(x):

$$\delta^2S=\int_{t_1}^{t_2}\left(m\,\dot\eta^2-V''(x)\,\eta^2\right)dt$$

> **Complemento** (a apostila afirma a estacionariedade; o exemplo abaixo é meu).

Teste com η = sen(π(t − t₁)/Δt), onde Δt = t₂ − t₁:

$$\delta^2S=\frac{\Delta t}{2}\left(\frac{m\pi^2}{\Delta t^2}-V''\right)$$

| Partícula parada em | V'' | Ação |
|---|---|---|
| **máximo** de V | negativa | δ²S > 0 sempre: **mínimo** local da ação |
| **mínimo** de V | positiva | deixa de ser mínimo quando Δt > π√(m/V'') |

Ou seja, mínimo do potencial e mínimo da ação são perguntas diferentes.

## 9. Atrito de Rayleigh

> **Complemento** (aparece nos Hands-On; a apostila não desenvolve).

Forças dissipativas não saem de um potencial. Para atrito proporcional à velocidade, usa-se R = ½ γ q̇²:

$$\frac{d}{dt}\frac{\partial L}{\partial\dot q}-\frac{\partial L}{\partial q}+\frac{\partial R}{\partial\dot q}=0$$

Derivando H e usando a equação acima:

$$\frac{dH}{dt}=-\dot q\,\frac{\partial R}{\partial\dot q}=-\gamma\,\dot q^2=-2R\le0$$

A energia só diminui, e a taxa é o dobro de R.

## Teste-se

<details>
<summary>Deduza Euler–Lagrange em quatro frases.</summary>

Perturbo a trajetória por εη com as pontas fixas. Derivo a ação em ε no zero. Integro por partes e o termo de borda some. Como o resultado deve ser zero para todo η, o fator que multiplica η é zero.
</details>

<details>
<summary>Quando H = q̇ ∂L/∂q̇ − L é constante?</summary>

Quando L não depende explicitamente do tempo (e não há dissipação).
</details>

<details>
<summary>Por que a corda sem restrição dá uma resposta "absurda"?</summary>

Porque sem fixar o comprimento ela pode se esticar sem limite para baixar a energia. O erro está na pergunta.
</details>

## Para ir mais fundo

- Exercícios da apostila, em especial 6.2 (Lagrangiana relativística) e 11.2 (a montanha-russa mais rápida).
- Goldstein, *Classical Mechanics*; Landau e Lifshitz, *Mechanics*; Gelfand e Fomin, *Calculus of Variations*.
- E. Noether, *Invariante Variationsprobleme* (1918).
