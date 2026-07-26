# Junction 辺コンダクタンス b_s(k)：分数次数ベッセルによる closed-out

`branching-gate-Jk.md` で `J_k` を circulant k端子ネットワークへ還元し、残る唯一の
未知量は junction の辺コンダクタンス `b_s(k)`（＝channel admittance `lambda_m`）と
した。本稿はそれを **modified Helmholtz `Δg = g/(4T^2)` の分数次数ベッセル解**へ帰着
させ、構造を厳密に定め、定数を数値で closed-out する。

## 1. 厳密なモード解（分数次数ベッセル）

gate 近傍 `F = -Re(w^k) = -ρ^k cos(kφ)`。共形写像 `W = w^k`（`U = Re W = ρ^k cos kφ`,
`R = |W| = ρ^k`, `Ψ = kφ`）で

```
div(e^{-F/T} ∇h) = 0   ==>   div(e^{U/T} ∇_W h) = 0   ==>   Δ_W h + (1/T) ∂_U h = 0.
```

`g = e^{U/(2T)} h` とおくと `Δ_W g = g/(4T^2)`（2 次元 modified Helmholtz, 遮蔽長 `2T`）。
`W`-polar での変数分離、`Ψ`-周期 `2πk`（k-sheeted cover）から角周波数は `m/k`
（`e^{i(m/k)Ψ} = e^{imφ}`, `m ∈ Z`）、動径は order `|m|/k` の変形ベッセル。原点で正則な解は

```
┌─────────────────────────────────────────────────────────────────────────┐
│  h(ρ,φ) = e^{-U/(2T)} Σ_m  α_m  I_{|m|/k}( ρ^k / (2T) )  e^{i m φ},         │
│  U = ρ^k cos(kφ).                                                          │
└─────────────────────────────────────────────────────────────────────────┘
```

**分数次数 `m/k`** が k 回分岐の署名。`K_{|m|/k}` は原点で発散するので落ちる
（gate 中心で `h` 有限）。

### 数値検証（`python/branching_bessel_modes.py`）

FD committor（k=4, T=0.2）を各 `e^{imφ}` に射影し `ĝ_m(ρ)`、`I_{|m|/k}(ρ^k/2T)` との比を
確認。`ρ ∈ [0.35, 0.95]`（well/confinement が効く前の gate 領域）で比はほぼ一定：

| mode m | order m/k | 比の相対ばらつき | 係数 α_m |
|---|---|---|---|
| 0 | 0 | 0.0000 | +0.5000 |
| 1 | 1/4 | 0.005 | +0.1433 |
| 2 | 1/2 | 0.007 | +0.2041 |
| 3 | 3/4 | 0.005 | +0.1433 |
| 4 | 1 | 0.0000 | +0.5000 |

分数次数ベッセルのモード構造が **1% 以内**で実証（ρ>0.9 の微小ドリフトは confinement
補正）。`α_1 = α_3`（A↔B 反射対称）、整数次数モード（m=0,4）は係数 1/2。

## 2. channel admittance と connection problem

circulant 固有ベクトル `c_j = e^{2π i m j/k}` の channel に admittance `λ_m` を対応。
谷 `φ_j = 2πj/k` での境界値 `e^{imφ_j}` は Fourier mode `m + nk`（`n ∈ Z`）を区別しない
（aliasing）。よって channel m は **整数間隔の order `|m/k + n|`** を励起する：

```
h^{(m)} = e^{-U/(2T)} Σ_n β_n I_{|m/k + n|}( ρ^k/(2T) ) e^{i(m+nk)φ}.
```

ridge（`cos kφ = -1`）上での boundedness（`h` 有界）がこの無限和の係数を選び、
admittance `λ_m` が定まる。この connection problem は一般に初等的閉形式を持たない。

### k=2（Morse）は厳密

`k=2`：`F = -(x^2 - y^2)`、非退化鞍点。committor は誤差関数
`h = (1/2) erfc(x/√(2T))`（不安定座標 x のみの関数）で、split 公式より `J_2 = 1`（厳密）。

### 数値 closed-out（k=3..6）

effective resistance `R_0d = 1/cap(0,d)` から circulant 固有値 `λ_m` を逆算し、
scale-invariant admittance `a_m = λ_m/T` と `J_k = 1/R_01` を得る（T=0.12）：

| k | a_1 | a_2 | a_3 | J_k = cap/T |
|---|---|---|---|---|
| 3 | 1.672 | — | — | 0.836 |
| 4 | 1.383 | 1.942 | — | 0.808 |
| 5 | 1.159 | 1.861 | — | 0.797 |
| 6 | 0.990 | 1.704 | 1.963 | 0.792 |

（`a_m = a_{k-m}`、独立成分のみ表示。）辺コンダクタンス `b_s` は `a_m` の逆 DFT。
`J_k` は §1 のネットワーク公式へ代入して閉じる。

## 3. 到達点

- **厳密（新）**：committor の分数次数ベッセルモード解 `h = e^{-U/2T} Σ α_m I_{|m|/k}(ρ^k/2T) e^{imφ}`。
  数値で 1% 以内検証。共形 → modified Helmholtz → 分数ベッセルの連鎖。
- **厳密**：`J_2 = 1`（Morse, erfc）。
- **数値的に確定**：channel admittance `a_m(k)` と `J_k`（k=3..6）。ネットワーク公式で
  `J_k` が閉じる。
- **残る解析課題**：admittance `a_m`（= `b_s` の Fourier 変換）の connection problem
  （aliased order `|m/k+n|` + ridge boundedness）の初等的閉形式。第一撃と同じ `K_0`/`I_ν`
  ベッセル族が支配する。`k→∞` 漸近も未了。

## 4. 証拠の区分

- **厳密**：モード解の分数次数ベッセル構造（変数分離の帰結）；`J_2 = 1`；
  channel と circulant 固有値の対応。
- **数値的に確定**：モード構造の射影検証（<1%）；`a_m(k)`, `J_k`（k=3..6）。
- **未決**：`a_m` の閉形式（Bessel connection problem）、`k→∞` 漸近、`b_s` の初等表示。
