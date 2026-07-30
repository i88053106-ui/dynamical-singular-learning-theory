# Case III：non-splittable gate と PDE-algebraic invariant

theory.md Case III：不安定方向と特異核が局所的に分離できず、単一の非退化 quadratic
unstable mode を split off できない gate。本稿は具体モデルで、**capacity が Newton
（RLCT）data では決まらず、germ の係数に依存する PDE/解析的不変量である**ことを示す。

## 1. モデル

```
F(x,y) = -x^4 + a x^2 y^2 + y^4 + (x^2+y^2)^3
```

- 井戸：±x 軸上 `x0 = sqrt(2/3)`（`F(x,0) = -x^4 + x^6`、a 非依存）。深さ `-4/27`。
- gate：原点、`H = 0`。degree-4 homogeneous germ `-x^4 + a x^2 y^2 + y^4`（confinement は degree 6）。
- 2 valley（±x）、非 branching。

**Case I・II の外にある理由。**

1. **Hessian = 0、不安定方向が quartic（`-x^4`）** → 非退化 quadratic unstable mode が無く、
   split-gate 公式（`cap = e^{-H/T} I_K sqrt(μT/2π)`, μ 必要）が適用不能。
2. **germ が非調和**（`Δ(-x^4+x^2y^2+y^4) = -10x^2+14y^2 ≠ 0`）→ branching（Case II）の
   共形線形化 `W = w^k` が適用不能。
3. **`a ≠ 0` で非分離**（`x^2 y^2` 結合）。`a = 0` の leading germ `-x^4 + y^4` は分離的で、
   committor は 1 次元（x のみ）に還元される。

## 2. scaling

homogeneous degree-4、dimension 2 なので `cap_T ~ J(a) T`（`p = 1 + (d-2)/m = 1`）。
数値（a=1、`python/nonsplittable_gate.py`）：

| T | cap | J = cap/T | p |
|---|---|---|---|
| 0.060 | 0.042914 | 0.7152 | |
| 0.045 | 0.032637 | 0.7253 | |
| 0.030 | 0.022401 | 0.7467 | |
| 0.020 | 0.015430 | 0.7715 | p ≈ 0.930 → 1 |

（barrier が浅く confinement 補正 ~`sqrt(T)` が大きいため有限 T の p は 1 よりやや小。）

## 3. 主結果：J は germ の係数に依存する（Newton 不変でない）

germ `-x^4 + a x^2 y^2 + y^4` の **Newton 多面体は `{(4,0),(2,2),(0,4)}` の凸包で、
係数 a にも符号にも依存しない**。それにもかかわらず capacity 定数 `J(a)` は a とともに
変化する（fixed T、外挿なし、両 T で単調）：

| a | J (T=0.045) | J (T=0.030) |
|---|---|---|
| 0.0 | 0.7994 | 0.8281 |
| 0.5 | 0.7602 | 0.7850 |
| 1.0 | 0.7253 | 0.7467 |
| 1.5 | 0.6940 | 0.7126 |
| 2.0 | 0.6660 | 0.6821 |
| 3.0 | 0.6181 | 0.6304 |

`a = 0 → 3` で `J` は **約 24% 減少**（両 T で同じ単調傾向、figures/nonsplittable_gate.png）。

```
┌────────────────────────────────────────────────────────────────────┐
│  Newton polygon 固定、係数 a のみ変化 → J(a) が変化.                   │
│  ゆえに capacity 定数は germ の PDE/解析的不変量であり、               │
│  Newton (RLCT) data の関数ではない.  [Case III]                       │
└────────────────────────────────────────────────────────────────────┘
```

Case I（split gate）では capacity data `(κ,r)` は横断 Newton data `(λ⊥,m⊥)` で決まった
（係数非依存）。Case II（branching）も homogeneous 構造から決まった。Case III では
**係数（不安定 vs 安定の符号、結合の強さ）が capacity に効く**——これが「PDE–algebraic
capacity invariant が必要」という theory.md の主張の具体的実証である。

## 4. 分離性の対比（a=0 vs a>0）

`a = 0`：leading germ `-x^4 + y^4` は分離的。分離ポテンシャル `-U(x) + V(y)` では
committor `h(x) = C ∫^x e^{-U/T}` が x のみの関数として解 `div(e^{(U-V)/T}∇h)=0` を満たす
（1 次元還元）。`a > 0`：`x^2 y^2` 結合が分離を壊し、committor は本質的に 2 次元、capacity は
full 2D solve を要する。J(a) の a 依存はこの非分離性の定量的帰結。

## 5. 証拠の区分

- **厳密（構造）**：Hessian=0（split 不適用）；非調和（共形線形化不適用）；`a=0` の
  leading germ の分離性と 1D committor 還元；`cap ~ J(a) T`（scaling, `p=1`）。
- **数値的に確定**：`J(a)` の a 依存（fixed T、両 T で単調、24% 変化）＝ capacity が
  Newton 不変でないこと。
- **未決**：`J(a)` の閉形式（新しい PDE-algebraic invariant）；一般の non-splittable
  germ に対する分類。これが完全 singular capacity theory の最終フロンティア。
