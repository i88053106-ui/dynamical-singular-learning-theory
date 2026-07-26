# Branching-gate 定数 J_k の解析

第三撃で数値確定した `cap_T ~ J_k T`（degree-k homogeneous branching gate、`p=1`）の
定数 `J_k` を解析的に構造化する。結論を先に述べる：

- **`J_k` は circulant な k端子ネットワークの隣接端子間 effective conductance**である
  （数値でパラメータフリーに検証、下記 §3）。
- したがって 2 次元 committor 問題は、**有限対称ネットワーク + junction の辺
  コンダクタンス `b_s(k)`** に厳密還元される。`b_s(k)` の閉形式が残る唯一の未解決点。
- gate の共形線形化（`W = w^k`）は重みを `e^{Re W}` に変え、operator を
  `Δh + ∂_U h = 0`（定数係数）へ帰着させる（§2）。

## 0. scale-invariant 変分による J_k

gate 近傍 `F ≈ -r^k cos(kθ) = -Re[(x+iy)^k]`（`H = 0`）。`x = T^{1/k} ξ` で
`F/T = -Re[(ξ_1+iξ_2)^k]` はスケール不変。confinement は inner で `O(T)` に落ちる。
capacity `cap_T = inf_h T ∫ |∇h|^2 e^{-F/T} dA` は

```
┌────────────────────────────────────────────────────────────────┐
│  cap_T = J_k · T,                                                 │
│  J_k = inf { ∫_{R^2} |∇h|^2 e^{ρ^k cos(kφ)} dA :                  │
│              h → 1 in valley A, h → 0 in valley B }.              │
└────────────────────────────────────────────────────────────────┘
```

谷の奥で committor の flux 保存から `|∇h|^2 e^{-F} ~ e^{F} → 0`、`J_k` は有限
（第三撃で `p=1`, log なしを数値確定）。

**k=2（Morse）アンカー。** `Re[w^2] = x^2 - y^2` は非退化鞍点。split 公式より
`cap = I_K sqrt(μ T/2π)`、`I_K = sqrt(πT)`, `μ = 2` で `cap = T`。よって `J_2 = 1`（厳密）。

## 1. 深い谷 = 端子：k端子ネットワークへの還元

各谷は深い井戸で、内部では committor がほぼ一定値 `h_j` に平衡化する。したがって
junction（gate 近傍）は **k 個の端子を持つ線形導体**であり、端子電位 `h_j` と注入電流
`I_j` は Laplacian（conductance）行列 `G` で `I = G h` と関係する。k 回対称性から
`G` は **circulant**：`G_{jl} = -γ_{|l-j|}`（`γ_s`= 距離 s の谷間 conductance,
`γ_s = γ_{k-s}`）、行和 0。

capacity は端子 0,1 間（他端子は float, `I_j = 0`）の effective conductance：

```
J_k = 1 / R_{01},   R_{01} = (1/k) Σ_{m=1}^{k-1} 4 sin^2(π m/k) / λ_m,
```

`λ_m = Σ_{s=1}^{k-1} γ_s (1 - cos(2π m s/k))` は `G` の固有値。spectator well の
committor 値は同じネットワークの floating-node 電位である。

### k=4 の明示形（単一結合比 r = c/b）

`γ_1 = b`（隣接）, `γ_2 = c`（対角）, `r = c/b`：

```
spectator (B 隣接) :  h = (1+r)/(3+r),        r = (3h-1)/(1-h),
cap(隣接)          :  J_4 = 4 b (1+r)/(3+r),
cap(対角)          :  b (1+r),
cap(対角)/cap(隣接) = (3+r)/4.
```

最後の比は base conductance `b` に依存せず、**spectator 値 `h` だけから決まる**。

## 2. 共形線形化：W = w^k

`div(e^{-F} ∇h) = 0`, `e^{-F} = e^{Re(w^k)}`。等角写像 `W = w^k`（Dirichlet 積分は
2 次元で共形不変）で `∫ |∇_w h|^2 e^{Re w^k} dA_w = ∫ |∇_W h|^2 e^{Re W} dA_W`。
重みは `e^{U}`（`U = Re W`）に化け、Euler–Lagrange は

```
div(e^{U} ∇h) = 0   <=>   Δh + ∂_U h = 0   (定数係数),
```

`g = e^{U/2} h` とおくと `Δg = g/4`（modified Helmholtz, 遮蔽長 2）。2 次元
modified Helmholtz の Green 関数は `K_0(|·|/2)/(2π)` で、第一撃の横断積分に現れた
同じ `K_0` である。k 個の谷は branched k-sheeted cover 上の境界条件になる。
junction 結合 `γ_s`（＝辺コンダクタンス `b_s`）はこの線形問題のシート間応答として
定まる——が、その閉形式評価は未了。

## 3. 数値検証（circulant 還元）

`python/branching_network_reduction.py`。

**k=4 パラメータフリー test：** spectator 値 `h` から `r = (3h-1)/(1-h)` を求め、
`(3+r)/4` が独立測定の `cap(対角)/cap(隣接)` に一致するか。

| T | cap_adj | cap_opp | h_spec(B) | r | pred (3+r)/4 | meas | 差 |
|---|---|---|---|---|---|---|---|
| 0.16 | 0.12810 | 0.10980 | 0.4167 | 0.4285 | 0.8571 | 0.8571 | 2e-8 |
| 0.13 | 0.10478 | 0.08973 | 0.4162 | 0.4257 | 0.8564 | 0.8564 | 3e-8 |
| 0.10 | 0.08113 | 0.06943 | 0.4157 | 0.4229 | 0.8557 | 0.8557 | 2e-7 |

**7–8 桁一致**。circulant 4 端子ネットワークが gate junction を厳密再現。

**spectator 対称性（k=3..6）：** 反射対称 `h → 1-h` により相補ペアは和 1、対称軸上は
1/2。`|Σ - 1|` の残差は 2e-8 (k=4) 〜 5e-5 (k=6) と数値誤差レベルで成立。

## 4. 到達点と残る未解決点

- **確定（数値）**：`cap_T = J_k T`、`J_k` = circulant k端子ネットワークの隣接
  conductance。k=4 のパラメータフリー関係を 7–8 桁で検証。spectator 分岐＝floating
  node 電位。`J_2 = 1`（厳密）。
- **構造（解析）**：共形写像 `W = w^k` が重みを `e^{Re W}` に線形化、`Δg = g/4`
  （K_0 Green 関数）。junction は branched cover 上のシート間応答。
- **残る唯一の未解決点**：junction 辺コンダクタンス `b_s(k)`（＝`γ_s`）の germ
  `Re[w^k]` からの閉形式評価。これが求まれば `J_k` は完全に閉じる
  （`J_4 = 4b(1+r)/(3+r)` 等に代入）。`b_s` は §2 の modified Helmholtz シート間
  応答として原理的に定義される。

## 5. 証拠の区分

- **厳密**：`J_2 = 1`；scale-invariant 変分の定式化；共形線形化 `Δh + ∂_U h = 0`；
  circulant effective-conductance 公式（グラフ Laplacian の標準結果）。
- **数値的に確定**：circulant 還元（k=4 パラメータフリー test 7–8 桁、spectator 対称性）。
- **未決**：`b_s(k)` の閉形式、`J_k` の閉形式、`k → ∞` 漸近。
