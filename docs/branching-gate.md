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
精密化（§3.5）で `p → 1`（log 補正なし）を確認。

### 三つの発見

1. **split 公式は原理的に適用不能**。gate の Hessian が零ゆえ不安定曲率 `mu` が
   存在しない。公式 `cap = e^{-H/T} I_K sqrt(mu T/2pi)` に `mu → 0` を入れると
   `cap → 0` を予測するが、真の capacity は有限（~0.16）。**定量的な破綻**。
2. **branching 署名 `h_C = 1/2`（厳密）**。第三の井戸 C は committor 0.5 の
   separatrix 上に乗る（A↔B 対称性 `h → 1-h` と C の固定点性から厳密に 1/2）。
   1 次元 A→B チャネルなら C は無関係のはず。committor が 2 次元であることの直接証拠
   （figures/branching_gate_committor.png：separatrix が原点を通り C を巻き込む）。
3. **有効指数 p ≈ 0.94**（低温で 1 へ収束、§3.5）。leading homogeneous scaling
   （`r = T^{1/3} rho` で `F/T = -rho^3 cos3θ` がスケール不変）は `cap ~ T^1` を
   予測。谷の奥では committor の flux 保存から `|∇h|^2 e^{-F} ~ e^{F} → 0` で
   scale-invariant 積分 `J` が収束するため、`p = 1`（log なし）が理論値。

## 3.5 精密化：p → 1 の確認

`python/branching_gate_refine.py`。格子収束（T=0.2 で n=281,361,481 の cap/T は
0.81534, 0.81521, 0.81512、相対幅 0.0003）で数値の信頼性を確認後、低温へ拡張。

| T | cap(A,B) | cap/T | 局所 slope p |
|---|---|---|---|
| 0.30 | 2.3792e-1 | 0.79308 | 0.9289 |
| 0.25 | 2.0086e-1 | 0.80343 | 0.9348 |
| 0.20 | 1.6304e-1 | 0.81521 | 0.9447 |
| 0.16 | 1.3205e-1 | 0.82533 | 0.9549 |
| 0.13 | 1.0830e-1 | 0.83311 | 0.9645 |
| 0.10 | 8.4090e-2 | 0.84090 | — |

- **局所 slope が単調に 1 へ**（0.929 → 0.965、T 減少で増加）。log 補正があれば
  この T 域で `p ≈ 0.4` になるはずで、明確に否定される。
- **`cap/T` は `T^{1/3}` の一次式**：`cap/T = J + b T^{1/3}`、`J = 0.952`,
  `b = -0.235`（figures/branching_gate_exponent.png、点がほぼ完全に直線上）。
  この `T^{1/3}` 補正形自体が degree-3 homogeneous scaling の直接の証拠。

結論：**monkey saddle capacity `cap_T ~ J·T`（`p = 1`, log なし, `J ≈ 0.95`）**、
有限 T 補正は予測どおり `T^{1/3}`。degree-3 branching gate の capacity 漸近が
数値的に確定した。

## 4. 意味と未解決問題

branching gate では、capacity が potential germ の RLCT/Newton data だけでなく
**committor 幾何・gate topology（谷の本数と配置）**に依存する。これが theory.md が
Case II/III として掲げた完全 singular capacity theory の中心的難所である。

**残る課題。**
- monkey saddle（degree-3）は `p = 1`（数値確定, §3.5）。定数 `J ≈ 0.95` の
  scale-invariant committor 積分による解析的評価は未了。
- k 本の谷（`Re[(s+iv)^k]`, degree-k homogeneous）への一般化。次元 2 での scaling
  予測は `cap ~ T^{1+(2-2)/k} = T^1`（k 非依存の power）だが、定数 `J_k` の k 依存性、
  および log 補正の有無は未検証。
- non-splittable gate（Case III：不安定方向と kernel が分離不能）への拡張。

## 4.5 谷数 k 依存性

`python/branching_valley_number.py`。degree-k homogeneous・k 回対称の gate 族

```
F_k(x,y) = r^{2k}/(2k) - Re[(x+iy)^k] = r^{2k}/(2k) - r^k cos(kθ),
```

k 個の井戸（θ_j = 2π j/k, r0 = k^{1/k}, 深さ -k/2）、共通 gate は原点、`H = 0`。
`k = 2` は普通の Morse 鞍点（split 公式が使える境界）、`k ≥ 3` が branching gate。
隣接する 2 井戸間の cap(A,B) を `k = 2..6` で測定。

| k | p_global | slope(low T) | cap/T @ T=0.10 | spectator well の committor |
|---|---|---|---|---|
| 2 | 0.914 | 0.940 | 0.951 | —（Morse） |
| 3 | 0.950 | 0.966 | 0.841 | {0.500} |
| 4 | 0.965 | 0.976 | 0.812 | {0.416, 0.584} |
| 5 | 0.973 | 0.981 | 0.800 | {0.384, 0.500, 0.616} |
| 6 | 0.977 | 0.980 | 0.795 | {0.388, 0.490, 0.544, 0.633} |

**発見。**

1. **power は k 非依存**：局所 slope はすべて 1 へ収束（0.94 → 0.98、k が大きいほど
   速い）。scaling 予測 `1 + (d-2)/k = 1`（d=2）を全 k で確認。
2. **capacity 定数は k とともに単調減少**：cap/T @ T=0.10 は 0.951 → 0.795
   （extrapolation-free、curves は panel A で非交差）。谷が増えるほど隣接対の
   直接 conductance は小さくなる。
3. **k=2 は Morse/split 境界**：cap/T → ~1（split 公式の leading `cap = T`）。
   branching solver が split 理論に接続することの確認。
4. **spectator committor の分岐構造**：傍観する k-2 個の井戸は k 回対称に従う特定の
   分数値を取る。反射対称 `h → 1-h` により相補ペアは和が 1
   （例 k=4: 0.416+0.584=1.000）、対称軸上の井戸はちょうど 1/2（k=3, k=5）。
   committor が 2 次元 harmonic であることの指紋。

**注意（正直な区分）。** 絶対正規化は committor solver の ~6% 系統誤差（Morse 検証）を
含む。また cap/T の有限 T 補正指数は k 依存（k=2 は `O(T)`、k=3 は `T^{1/3}`）で
`T→0` 外挿はモデル依存。よって k 間比較は同一 T・同一 solver の cap/T（外挿なし）で
行い、絶対 `J_k` の精密値は主張しない。robust なのは「p_k=1」「k で単調減少」
「spectator 分岐構造」の 3 点。

## 5. 証拠の区分

- **厳密/既知**：monkey saddle の臨界構造、`H = 0`、`h_C = 1/2` の対称性論法、
  Hessian = 0 による split 公式の不適用。
- **数値的に確定**：solver 検証（Morse ~6%、格子収束 0.03%）、`h_C = 0.5000`、
  cap の指数 `p → 1`（局所 slope 0.929→0.965、`cap/T` が `T^{1/3}` 一次で `J≈0.95` へ）。
- **予想/未決**：定数 `J` の解析評価、谷数 k への一般化、branching-gate capacity の
  分類不変量。
