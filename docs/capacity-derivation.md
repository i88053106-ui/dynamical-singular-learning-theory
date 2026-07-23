# Split-gate capacity の導出

第一標的モデルにおける transition capacity の漸近を、検証済みの横断特異積分
`I(T)` に還元して導く。ここでの導出は **formal（sharp な定数まで込みの漸近）**
であり、rigorous な下界（Newton-adapted lower bound, theory.md §9）は別途の
定理課題として残る。

## 1. 設定と規約

可逆 overdamped Langevin と非正規化 Dirichlet form（theory.md §3）：

```
E_T(h,h) = T ∫ |∇h(w)|^2 e^{-F(w)/T} dw,
cap_T(A,B) = inf { E_T(h,h) : h|_A = 1, h|_B = 0 }.
```

第一標的 potential：

```
F(x,y,z) = (x^2 - 1)^2/4 + x^2 (y^2 + z^2) + y^6 + y^2 z^2 + z^6,
K(y,z)   = y^6 + y^2 z^2 + z^6.
```

wells は `(±1, 0, 0)`（Morse, F=0, Hess=diag(2,2,2), det=8）、
gate は原点、communication height `H = 1/4`。

**厳密な split normal form（恒等式）。** `X = x·sqrt(1 - 2(y^2+z^2) - x^2/2)` と
おくと、原点近傍で

```
F - 1/4 = -X^2/2 + K(y,z)            (exact identity)
```

が成り立つ。実際 `-X^2/2 = -x^2/2 + x^2(y^2+z^2) + x^4/4` であり、
`F - 1/4 = -x^2/2 + x^2(y^2+z^2) + x^4/4 + K` に一致する。よって不安定曲率は
`μ = 1`、横断幾何は Newton-degenerate な `K`。

## 2. Split-gate capacity 公式

ゲートが厳密に `F - H = -(μ/2)s^2 + K(v)`（`s` 不安定 1 次元、`v` 横断）へ
分離するとき、equilibrium potential は主要オーダーで `s` のみに依存する。
このとき Dirichlet energy は

```
E_T(h,h) = T e^{-H/T} [∫ e^{-K(v)/T} dv] [∫ |h'(s)|^2 e^{(μ/2)s^2/T} ds]
```

と因子化し、1 次元変分 `inf_h ∫ |h'|^2 e^{(μ/2)s^2/T} ds = [∫ e^{-(μ/2)s^2/T} ds]^{-1}
= sqrt(μ/(2πT))` を用いると

```
┌─────────────────────────────────────────────────────────┐
│  cap_T(A,B) = e^{-H/T} · I_K(T) · sqrt(μ T / (2π)),       │
│  I_K(T) = ∫ e^{-K(v)/T} dv.                                │
└─────────────────────────────────────────────────────────┘
```

### 古典極限とのクロスチェック（数値で機械精度一致を確認済み）

- **1 次元 Kramers。** 横断次元 0（`I_K = 1`）で
  `cap_T = e^{-H/T} sqrt(μT/2π)`。`E[τ] = Z_well/cap` は古典 Kramers
  `2π/(sqrt(F''(m)) sqrt(μ)) e^{ΔH/T}` を厳密再現。
- **多次元 Morse。** `K = Σ_{i≥2} (μ_i/2) v_i^2` のとき
  `I_K = ∏ sqrt(2πT/μ_i)`。`E[τ] = Z_well/cap` は標準 Eyring–Kramers
  `(2π/μ) sqrt(|det Hess F(σ)| / det Hess F(m)) e^{ΔH/T}` を厳密再現。

したがって公式 `cap_T = e^{-H/T} I_K sqrt(μT/2π)` は Eyring–Kramers 則の
（Newton-degenerate 横断への）正しい一般化である。

## 3. 本モデルへの適用

`H = 1/4`, `μ = 1`, `I_K(T) = I(T)`（検証済み横断積分）。
`I(T) = (sqrt(π)/3) sqrt(T) [log(1/T) + γ + 6 log 2] + O(T^{3/2} log(1/T))`
を代入し `sqrt(π)/sqrt(2π) = 1/sqrt(2)` を使うと：

```
┌───────────────────────────────────────────────────────────────────┐
│  cap_T(A,B) = (1/(3 sqrt2)) T [log(1/T) + γ + 6 log 2] e^{-1/(4T)}  │
│             + O(T^2 log(1/T)) e^{-1/(4T)}.                          │
└───────────────────────────────────────────────────────────────────┘
```

leading:  `cap_T ~ (1/(3 sqrt2)) T log(1/T) e^{-1/(4T)}`（theory.md §8/§19 の予想 A と一致）。

数値確認（e^{-1/4T} を除いた比 cap/refined）：T=1e-2 で 1.00503、T=1e-3 で 1.00073
→ `O(T)` で 1 へ収束。

## 4. Exit-time への帰結

`Z_A = ∫_A e^{-F/T} ~ (2πT)^{3/2}/sqrt(8)`（Morse well, det=8, F(m)=0）。
potential-theoretic 公式 `E[τ_{A→B}] ~ Z_A / cap_T(A,B)` より：

```
┌────────────────────────────────────────────────────────────────┐
│  E[τ_{-→+}] ~ 3 sqrt2 · π^{3/2} · T^{1/2} / log(1/T) · e^{1/(4T)} │
└────────────────────────────────────────────────────────────────┘
```

対数表示：

```
log E[τ] = 1/(4T) + (1/2) log T - log log(1/T) + O(1).
```

すなわち negative logarithmic prefactor exponent `(log(1/T))^{-1}` が
Newton-degenerate transverse gate の観測署名（theory.md §8）。

## 5. 証拠の区分と残る課題

- **厳密（恒等式・既知理論）**：split normal form の恒等式；公式が 1D Kramers /
  多次元 Eyring–Kramers を再現すること。
- **formal（sharp 漸近, 数値で強く支持）**：本モデルの capacity 定数 `1/(3 sqrt2)`
  と refined 補正 `γ + 6 log 2`。upper bound は test 関数で容易。
- **残る定理課題（theory.md §9）**：matching lower bound を Newton-adapted
  sublevel section `D_T = {K < M T log(1/T)}` 上で構成し、
  fiberwise reduction の Jacobian / metric 補正が `1 + o(1)` であることを示す。
- **未検証（次の数値課題）**：3 次元 Langevin の平均脱出時間 `E[τ]` を直接
  シミュレーションし、上記 exit-time 則（特に `(log(1/T))^{-1}` 署名）を確認する。
