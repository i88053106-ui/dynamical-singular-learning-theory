# 多価 non-splittable germ：係数依存の branching admittance

Case II（branching, 調和 `Re w^k`）は普遍 admittance `a_m = 2 sin(πm/k)`（係数非依存）、
Case III（non-splittable）は係数依存の capacity を与えた。両者を合流させ、**branching
network を保ちつつ admittance が係数依存になる**（＝ PDE-algebraic）ことを示す。

## 1. モデル：4 谷 branching gate

同次 degree-4 の族

```
P_a(x,y) = -x^4 + a x^2 y^2 - y^4      (a > 2)
```

は ±x, ±y の **4 谷**（対角線 x=±y は `(a-2)x^4 > 0` の稜線）、原点にゲート、4 回対称
（circulant 4 端子 junction ＝ Case II 構造）。**調和になるのは a=6 のみ**：

```
P_6 = -(x^4 - 6 x^2 y^2 + y^4) = -Re(w^4),   Δ P_a = (2a-12)(x^2+y^2)  ⇒  調和 ⟺ a=6.
```

- **a=6（調和）**：k=4 branching gate `-Re w^4`。scale-invariant admittance は普遍
  Steklov 値 `a_m = 2 sin(πm/4)`（`a_1 = √2 = 1.4142`, `a_2 = 2`）。
- **a≠6（非調和）**：admittance `a_m(a)` は **a に依存**（Case III）。

## 2. 数値（`python/branching_nonsplittable.py`）

`P_a` は同次なので scale-invariant `a_m` を純 germ・T=1・4 reservoir（谷の奥）で直接計算。

**a=6 の収束**（reservoir を深く → `a_1 → 2 sin(π/4)`）：

| L, n | a_1 |
|---|---|
| 1.3, 421 | 1.5993 |
| 1.6, 481 | 1.4535 |
| 1.9, 541 | 1.4173 |

→ `2 sin(π/4) = 1.4142`（0.2% 一致）。

**admittance vs 結合 a**（L=1.9, n=541）：

| a | a_1 | a_2 |
|---|---|---|
| 3 | 2.4659 | 4.0064 |
| 4 | 1.8623 | 2.8348 |
| 5 | 1.5866 | 2.3156 |
| **6** | **1.4173** | **2.0059** |
| 7 | 1.2985 | 1.7946 |
| 8 | 1.2087 | 1.6385 |

`a_1(a)`, `a_2(a)` はそれぞれ Steklov 線 `2 sin(π/4)`, `2 sin(π/2)` を **ちょうど a=6 で
横切り**、離れると偏差する（figures/branching_nonsplittable.png）。

## 3. 合流の言明

```
┌────────────────────────────────────────────────────────────────────┐
│  多価 non-splittable gate では:                                       │
│   ・branching network 構造（circulant k 端子）は保たれる（Case II）,  │
│   ・admittance a_m(a) は germ 係数に依存する PDE-algebraic 量（Case III),│
│   ・普遍 Steklov 値 a_m = 2 sin(πm/k) は調和点でのみ実現される.        │
└────────────────────────────────────────────────────────────────────┘
```

Case II の普遍性（`2 sin(πm/k)`）は germ が **調和（Re w^k）** であることに由来し、
非調和な多価 germ では admittance が係数依存になる。branching と non-splittable は
「network は共通・admittance の普遍性は調和性に依存」という形で統一される。

## 4. 証拠の区分

- **厳密/既知**：`P_a` の 4 谷・4 回対称・調和条件 `a=6`；調和点の Steklov 値
  `2 sin(πm/4)`（Case II）。
- **数値的に確定**：`a=6` で `a_m → 2 sin(πm/4)`（収束）；`a≠6` で `a_m(a)` の係数依存
  （a_1: 2.47→1.21, a_2: 4.01→1.64 over a=3..8）。
- **未決**：`a_m(a)` の閉形式（Case III 同様に非初等の見込み）；一般 k・一般非調和多価
  germ の admittance 分類；調和点まわりの摂動係数 `da_m/da|_{a=6}` の解析評価。
