# Split-gate capacity の一般分類（第二撃）

theory.md 第二撃・予想 B の内容を、capacity 公式（docs/capacity-derivation.md）と
横断特異 Laplace 幾何の接続として定式化し、germ ライブラリで数値検証する。

## 1. 分類定理候補

**設定。** gate が局所的に **index-one split form**

```
F - H = -(mu/2) s^2 + K(v) + R(s,v),    mu > 0,
```

（`s`：非退化な 1 次元不安定方向、`v ∈ R^{d-1}`：横断方向、`R` は主部より高次）へ
還元されるとする。横断 germ `K` は Newton-nondegenerate な非負解析関数で、その
横断 Laplace 積分が

```
I_K(T) = int e^{-K(v)/T} dv ~ C_K T^{lambda_perp} (log 1/T)^{m_perp - 1}
```

を満たすとする。

**主張（split-gate classification）。** capacity 公式
`cap_T = e^{-H/T} I_K(T) sqrt(mu T / 2pi)` に代入すると

```
┌────────────────────────────────────────────────────────────────────┐
│  cap_T(A,B) ~ e^{-H/T} · C_K sqrt(mu/2pi) · T^{lambda_perp + 1/2}     │
│                                       · (log 1/T)^{m_perp - 1},       │
│                                                                      │
│  すなわち  (kappa, r) = (lambda_perp + 1/2, m_perp).                   │
└────────────────────────────────────────────────────────────────────┘
```

capacity の power exponent は **横断 RLCT + 1/2**（`+1/2` は不安定方向の
Gaussian 寄与）、log multiplicity は **横断多重度そのもの**。

**2 変数横断の Newton 読み取り。** `d` を Newton 距離（対角線が Newton 多面体
境界に当たる点の `(d,d)`）とすると

```
lambda_perp = 1/d,
m_perp = 2  if 対角線が頂点（0 次元面）に当たる,
         1  if 辺（1 次元面）の内部に当たる.
```

一般に `m_perp = (d_perp) - k`（`k` = 対角線が当たる面の次元、`d_perp` = 横断次元）。

## 2. germ ライブラリによる数値検証

`python/split_gate_classification.py`。各 germ の monomial support から Newton
pole data を計算し、独立に `I_K(T)` を 2 次元 log-map quadrature で数値評価。
`r(T) = I_K(T) / [T^{lambda_perp} (log 1/T)^{m_perp-1}]` が定数へ収束すれば予測を追認。

| germ K(y,z) | Newton (λ⊥, m⊥) | 予測 (κ, r) | r(T): T=1e-2, 1e-3, 1e-4 |
|---|---|---|---|
| y² + z² | (1.000, 1) | (1.500, 1) | 3.14159, 3.14159, 3.14159 |
| y⁴ + z² | (0.750, 1) | (1.250, 1) | 3.21311, 3.21311, 3.21311 |
| y⁴ + z⁴ | (0.500, 1) | (1.000, 1) | 3.28626, 3.28626, 3.28626 |
| y⁶ + z⁶ | (0.333, 1) | (0.833, 1) | 3.44265, 3.44265, 3.44265 |
| y⁴ + y²z² + z⁴ | (0.500, 1) | (1.000, 1) | 2.98791, 2.98791, 2.98791 |
| y⁶ + y²z² + z⁶ | (0.500, 2) | (1.000, 2) | 1.20446, 0.99662, 0.89471 |

- `m=1` の germ は `r(T)` が全 T で**厳密に一定**（log 因子なし）。
  `y^2+z^2` は `r = pi`（`I_K = pi T` の厳密値）を機械精度で再現。
- `y^6+y^2z^2+z^6`（`m=2`）は `r(T)` が `sqrt(pi)/3 = 0.5908` へ収束
  （ドリフトは `gamma + 6 log2` 補正、`r = (sqrt(pi)/3)(1 + (gamma+6log2)/log(1/T))`）。

**決定的対比。** `y^4 + y^2 z^2 + z^4` と `y^6 + y^2 z^2 + z^6` は見かけ上どちらも
「対角の平方項 + 中央結合項」だが：

- quartic：中央項 `(2,2)` は辺 `(4,0)-(0,4)` 上に**共線**（Newton 図の頂点でない）
  → `m=1`、**log なし**。`I_K/T^{1/2}` は一定（2.988）。
- sextic：中央項 `(2,2)` は Newton 図の**真の頂点** → `m=2`、**log あり**。
  `I_K/T^{1/2}` は `log(1/T)` として増大（5.54 → 6.88 → 8.24）。

log multiplicity は Hessian やナイーブな次数勘定では決まらず、**Newton 頂点条件**が
支配する。これが split-gate capacity 分類の要点である。

## 3. metastable singularity gap への帰結

well 側 pole data `(lambda_A, m_A)` と合わせ、transition-edge data は

```
delta_{A->B} = lambda_A - kappa = lambda_A - (lambda_perp + 1/2),
q_{A->B}     = m_A - r        = m_A - m_perp.
```

第一標的（Morse well: `lambda_A = 3/2, m_A = 1`; sextic gate: `kappa = 1, r = 2`）では

```
delta = 3/2 - 1 = 1/2,   q = 1 - 2 = -1,
```

すなわち `E[tau] ~ T^{1/2} (log 1/T)^{-1} e^{Delta H/T}`（負の log exponent、既検証）。

## 4. 高次予測（未数値検証）

横断次元 `d_perp = 3`（ambient `d = 4`）で対角線が Newton 頂点に当たる germ、
例 `y^6 + z^6 + w^6 + (mixed vertex terms)` では `m_perp = 3` となり、capacity は
`(log 1/T)^2` を持つ。分類式 `(kappa, r) = (lambda_perp + 1/2, m_perp)` の直接予測。

## 5. 証拠の区分と残る課題

- **厳密/既知**：capacity 公式の古典極限一致（docs/capacity-derivation.md）；
  Newton-nondegenerate germ の `I_K` pole data（Arnold–Varchenko–Watanabe）。
- **数値的に確定**：6 germs の `(lambda_perp, m_perp)` と `(kappa, r)` の対応；
  vertex/edge による log の有無（quartic vs sextic の決定的対比）。
- **残る定理課題**：
  1. remainder `R(s,v)` を含む split form での capacity 上下界の sharp 化
     （theory.md §9 の Newton-adapted lower bound）。
  2. 多次元不安定（index ≥ 2）や退化不安定方向への拡張。

## 6. 次のフロンティア：branching gate（第三撃）

split classification は **unstable mode を横断特異幾何から分離できる**ことに依存する
（theory.md Case I）。これが壊れるのが **branching gate**（Case II）である。

典型例は monkey saddle

```
F ~ H + Re[(s + i v)^3] = H + s^3 - 3 s v^2 + ...
```

原点で **3 つの下り谷**が 120 度間隔で出会う。committor は単純な `0/1` boundary
layer では決まらず、複数 valley 間の **harmonic splitting**（真に 2 次元の committor
幾何）が必要。したがって capacity は「1 次元 conductance × 横断質量」へ因子化せず、

```
cap_T != e^{-H/T} I_K sqrt(mu T/2pi)
```

となる可能性が高い。ここで初めて potential germ だけでなく **committor germ /
gate topology** が capacity を決める。これが完全な singular capacity theory の
中心的未解決問題である。

第一段階の課題：monkey-saddle gate で committor PDE を有限差分で解き、capacity が
ナイーブ split 公式からどれだけ乖離するかを数値的に測定する。
