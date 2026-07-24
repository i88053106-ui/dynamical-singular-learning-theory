# Branching gate（第三撃）：split 分解が壊れる世界

split-gate classification（docs/split-gate-classification.md）は、gate を
`F - H = -(mu/2) s^2 + K(v)` へ **分離できる**ことに依存する。この前提が壊れるのが
branching gate（theory.md Case II）である。ここでは committor が単純な 1 次元境界層で
決まらず、potential germ だけでは capacity が定まらない。

## 1. モデル：monkey saddle 三井戸

```
F(x,y) = (1/6)(x^2 + y^2)^3 - (x^3 - 3 x y^2)
       = (1/6) r^6 - r^3 cos(3θ).
```

- 3 つの井戸：θ = 0, 120, 240 度、`r0 = 3^{1/3} ≈ 1.442`、`F = -3/2`。
- gate：原点。局所的に `F ≈ -r^3 cos(3θ)`（degree-3 homogeneous, **Hessian = 0**）。
  3 つの下り谷が 120 度間隔で出会う monkey saddle。
- 任意の井戸対を隔てる最低鞍点はすべて原点（`θ = 60,180,300` 稜線上で
  `F = (1/6) r^6 + r^3` は `r=0` 最小）。communication height `H = F(0) = 0`。

したがって A→B の遷移は、第三の谷 C が同じ原点で開いている中を通る。committor は
本質的に 2 次元。

## 2. 数値 committor solver

`python/branching_gate_committor.py`。自己随伴形 `div(e^{-F/T} grad h) = 0` を
finite-volume 離散化（対称正定値、upwind 不要）。`h = 1` on A、`h = 0` on B、外周は
no-flux Neumann。capacity は離散最小化子の Dirichlet energy

```
cap_T(A,B) = T * sum_edges a_edge (h_i - h_j)^2,   a = e^{-F/T}.
```

### solver 検証（Morse 二重井戸）

`F = (x^2-1)^2/4 + y^2/2`（omega=1）。split 公式は `cap_T = (T/omega) e^{-1/(4T)}`。

| T | cap 数値 | cap 公式 | 相対差 |
|---|---|---|---|
| 0.10 | 7.624e-3 | 8.209e-3 | 0.071 |
| 0.08 | 3.287e-3 | 3.515e-3 | 0.065 |
| 0.06 | 8.811e-4 | 9.302e-4 | 0.053 |

数値は公式を ~5–7% で再現（有限格子・有限井戸半径による系統オフセット、T 減少で
相対差も減少 → 指数スケーリングは正しい）。solver は妥当。

## 3. branching gate の結果

| T | cap(A,B) | h_C（予測 0.5） |
|---|---|---|
| 0.30 | 2.379e-1 | 0.5000 |
| 0.25 | 2.009e-1 | 0.5000 |
| 0.20 | 1.631e-1 | 0.5000 |
| 0.15 | 1.242e-1 | 0.5000 |

有効指数：`cap ~ C T^p`, **p ≈ 0.94**（`H = 0` なので指数因子 `e^{-H/T} = 1`）。

### 三つの発見

1. **split 公式は原理的に適用不能**。gate の Hessian が零ゆえ不安定曲率 `mu` が
   存在しない。公式 `cap = e^{-H/T} I_K sqrt(mu T/2pi)` に `mu → 0` を入れると
   `cap → 0` を予測するが、真の capacity は有限（~0.16）。**定量的な破綻**。
2. **branching 署名 `h_C = 1/2`（厳密）**。第三の井戸 C は committor 0.5 の
   separatrix 上に乗る（A↔B 対称性 `h → 1-h` と C の固定点性から厳密に 1/2）。
   1 次元 A→B チャネルなら C は無関係のはず。committor が 2 次元であることの直接証拠
   （figures/branching_gate_committor.png：separatrix が原点を通り C を巻き込む）。
3. **有効指数 p ≈ 0.94**。leading homogeneous scaling（`r = T^{1/3} rho` で
   `F/T = -rho^3 cos3θ` がスケール不変）は `cap ~ T^1` を示唆するが、下り谷方向で
   `e^{-F/T}` が増大し `r^6` confinement で cut off されるため、有限 T では
   `p < 1` の crossover。真の漸近不変量の決定が **開かれた問題**。

## 4. 意味と未解決問題

branching gate では、capacity が potential germ の RLCT/Newton data だけでなく
**committor 幾何・gate topology（谷の本数と配置）**に依存する。これが theory.md が
Case II/III として掲げた完全 singular capacity theory の中心的難所である。

**残る課題。**
- monkey saddle capacity の sharp 漸近不変量の同定（`p = 1` か、log 補正や
  他の指数か）。committor の scale-invariant germ の解析が必要。
- k 本の谷（`Re[(s+iv)^k]`）への一般化と、谷数 k に対する capacity 依存性。
- non-splittable gate（Case III：不安定方向と kernel が分離不能）への拡張。

## 5. 証拠の区分

- **厳密/既知**：monkey saddle の臨界構造、`H = 0`、`h_C = 1/2` の対称性論法、
  Hessian = 0 による split 公式の不適用。
- **数値的に確定**：solver 検証（Morse ~6%）、`h_C = 0.5000`、cap 数値と `p ≈ 0.94`。
- **予想/未決**：monkey saddle capacity の真の漸近指数（leading argument は `T^1`）、
  branching-gate capacity の分類不変量。
