# Aula de 16/11: cálculo das variações e mínima ação

[Índice](../README.md) · [Teoria](../teoria/02_calculo_variacional.md) · [Hands-On](../atividades/topico02_como_explicar.md) · [Próxima aula](2026-11-18_pendulo.md)

| Professor | Tópico | Material do email |
|---|---|---|
| Reinaldo Rosa | 2 | slides (17 páginas) e apostila *The Variational Principle* |

## Em uma frase

Em vez de escrever forças vetor por vetor, escrevemos **um único número para cada trajetória** (a ação) e pedimos que ele seja estacionário. As equações de movimento saem daí.

## O fio da aula

**1. Ponte com o Tópico 1.** Grandezas físicas são tensores de várias ordens.

| Ordem | Exemplo | Relação |
|---|---|---|
| 0 | energia potencial | V |
| 1 | força | F = m a |
| 2 | tensão | σ |
| 3 | piezoelétrico | d |
| 4 | rigidez elástica | σ = C ε |

**2. Do tensor ao escalar.** O slide resume como `3² → 3¹ → 3⁰` (em 3D, 9, 3 e 1 componentes): forças → momento → trabalho ≈ ação. A formulação variacional **não elimina** a estrutura tensorial; ela a codifica em um funcional escalar.

**3. Quatro formulações da mesma mecânica.**

| Ano | Autor | Objeto central |
|---|---|---|
| 1687 | Newton | F = m a |
| 1788 | Lagrange | L = T − V |
| 1834–35 | Hamilton | H = Σ pᵢ q̇ᵢ − L |
| anos 1830 | Jacobi | ∂S/∂t + H = 0 |
| 1918 | Noether | simetria ⇒ conservação |

**4. A matemática.** Mínimo de função (f'(x) = 0), variação de uma função, funcional, e a condição δS = 0 que leva a Euler–Lagrange.

**5. Os três Hands-On** (páginas 14 a 16 dos slides, iguais às páginas 32 a 34 da apostila).

## Ideias para guardar

- **Ação não é trabalho.** Ela tem dimensão de energia × tempo (J·s).
- **δ não é d.** `dy` é a mudança ao andar sobre a mesma curva; `δy` é a diferença entre duas curvas no mesmo x.
- **Estacionário não é o mesmo que mínimo.** Pode ser mínimo, máximo ou sela.
- **Mínimo de V** indica equilíbrio estável; **máximo de V**, instável.

| Simetria da ação | Quantidade conservada |
|---|---|
| translação no tempo | energia |
| translação no espaço | momento linear |
| rotação | momento angular |

## Cuidados com a notação do slide

> Observações de estudo, a confirmar com o professor.

- Na página 13 aparece `y(x) = y(x) + δy(x)`. Leia como "curva nova = curva antiga + δy". Com o mesmo símbolo dos dois lados, a conta só faz sentido pelas cores.
- Na tabela da página 2, o piezoelétrico aparece como ε_ij = d_ijk E_k, mas a densidade de energia usa a ordem d_ijk E_i ε_jk. Confira em qual índice o campo elétrico entra.

## Teste-se

<details>
<summary>O que é um funcional? Dê um exemplo.</summary>

É uma regra que recebe uma função inteira e devolve um número. Exemplo: o comprimento de uma curva, ou a ação S[q].
</details>

<details>
<summary>Por que se exige que a variação seja zero nas pontas?</summary>

Porque comparamos trajetórias que saem e chegam nos mesmos pontos. Isso zera o termo de borda na integração por partes.
</details>

<details>
<summary>Onde está o mínimo de f(x) = x² + 2x + 5 (exemplo do slide)?</summary>

f'(x) = 2x + 2 = 0 dá x = −1. De fato f(x) = (x + 1)² + 4.
</details>

## Minhas anotações

<!-- Escreva aqui o que lembrar da aula. -->
