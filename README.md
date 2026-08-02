# 動的特異学習理論 / Dynamical Singular Learning Theory (DSLT)

> **井戸の特異幾何は質量を分類する。ゲートの特異幾何は capacity（伝導度）を分類する。
> 準安定時間は両者の比として現れる。**

渡辺澄夫の特異学習理論（SLT）が静的に分類する特異 Laplace 幾何を、可逆 Langevin 系の
**質量–capacity 構造**へ接続し、学習ダイナミクスの準安定遷移時間を記述する研究プログラム。

可逆 overdamped Langevin `dW = -∇F dt + √(2T) dB` に対し、非正規化 Dirichlet form
`E_T(h,h) = T ∫ |∇h|² e^{-F/T}` の下で

```
E[τ_{A→B}] ~ well mass Z_A(T) / transition capacity cap_T(A,B).
```

中心的問い：**SLT が singular Laplace integral を RLCT pole data で分類するように、
singular transition capacity を分類する幾何学的不変量を構成できるか。**

---

## 現在の到達点：三つの標的

| 標的 | 内容 | 状態 | 主な文書 |
|---|---|---|---|
| **第一撃** | Newton-degenerate transverse gate | ✅ 解析導出＋数値検証（Julia/Python/Langevin） | [`docs/capacity-derivation.md`](docs/capacity-derivation.md) |
| **第二撃** | split gate の一般分類 `(κ,r)=(λ⊥+½, m⊥)` | ✅ 6 germ で検証 | [`docs/split-gate-classification.md`](docs/split-gate-classification.md) |
| **第三撃** | branching gate（monkey saddle） | ✅ split 破綻を実証・`p=1` 確定・谷数 k 依存を解明 | [`docs/branching-gate.md`](docs/branching-gate.md) |

理論全体の設計は [`docs/theory.md`](docs/theory.md)。

---

## 第一撃：Newton-degenerate transverse gate

第一標的 potential `F = (x²-1)²/4 + x²(y²+z²) + y⁶+y²z²+z⁶`、gate は原点、`H = 1/4`。
横断 germ `K = y⁶+y²z²+z⁶`。

**横断特異積分**（厳密な 1 次元 Bessel 表示 `I = (4/3)√T ∫ e^{-x²} K₀(2√T x³) dx` から）：

```
I(T) = ∫∫ e^{-K/T} dy dz = (√π/3) √T [ log(1/T) + γ + 6 log 2 ] + O(T^{3/2} log(1/T)).
```

独立 2 次元積分（HCubature/scipy）と 1 次元表示が相対差 ~1e-10 で一致（RLCT `λ=1/2`,
multiplicity `m=2`）。

**capacity と exit-time**（split-gate 公式 `cap_T = e^{-H/T} I_K √(μT/2π)`, 古典
Eyring–Kramers を機械精度で再現）：

```
cap_T ~ (1/(3√2)) T [log(1/T)+γ+6log2] e^{-1/(4T)},
E[τ]  ~ 3√2 π^{3/2} · √T / log(1/T) · e^{1/(4T)}.
```

負の logarithmic prefactor `(log(1/T))^{-1}` が観測署名。**Langevin シミュレーション**で
7 温度にわたり確認：`C_eff = log E[τ] - 1/(4T) - ½log T + log log(1/T)` が
**2.055 ± 0.025** と安定（`-log log(1/T)` 項の存在）。`1/log(1/T)` 因子を入れると
reduced 量は平坦、入れないと単調ドリフト（識別テスト）。

## 第二撃：split gate の一般分類

index-one split gate `F - H = -(μ/2)s² + K(v) + …` で、横断 germ の Laplace pole data
`(λ⊥, m⊥)` が capacity を決める：

```
(κ, r) = (λ⊥ + ½, m⊥).
```

`λ⊥ = 1/d`（Newton 距離）、`m⊥ = 2`（対角線が Newton 頂点）or `1`（辺）。6 germ で検証。
決定的対比：`y⁴+y²z²+z⁴`（中央項が辺上に**共線** → `m=1`, log なし）vs
`y⁶+y²z²+z⁶`（中央項が**頂点** → `m=2`, log あり）。log multiplicity は Hessian でなく
**Newton 頂点条件**が支配する。

## 第三撃：branching gate

split 分解が壊れる gate。monkey saddle `F = (1/6)r⁶ - (x³-3xy²)`（3 井戸、原点で 3 谷が
出会う、Hessian = 0）：

- **split 公式は原理的に不適用**（不安定曲率 μ が存在しない。μ→0 は cap=0 を予測するが
  真値は有限）。
- **branching 署名 `h_C = 1/2`**（第三の井戸が committor 0.5 の separatrix 上）。
  committor が 2 次元 harmonic である直接証拠。
- **capacity 指数 `p = 1`（log なし）を確定**：局所 slope 0.929→0.965 が 1 へ、
  `cap/T` が `T^{1/3}` の直線で `J ≈ 0.95` へ外挿。degree-3 homogeneous scaling の予測通り。
- **谷数 k 依存性**（`Re[(x+iy)^k]`, k=2..6）：power は k 非依存（全 k で `p→1`）、
  capacity 定数は k とともに単調減少、spectator well は k 回対称に従う分数 committor
  （相補ペアの和=1、対称軸上は½）。`k=2` は Morse/split 境界。

---

## リポジトリ構成

```
docs/        理論と各撃の導出
  theory.md                     全体設計（三層構造・予想）
  capacity-derivation.md        第一撃：split-gate capacity 公式
  split-gate-classification.md  第二撃：横断 Newton pole data による分類
  branching-gate.md             第三撃：monkey saddle・谷数 k 依存
  branching-gate-Jk.md          第三撃：定数 J_k の circulant ネットワーク還元
  branching-gate-bs.md          第三撃：b_s の分数次数ベッセルモード解
  branching-gate-largek.md      第三撃：admittance の k→∞ 漸近 a_m ~ 2πm/k
  branching-gate-slope.md       第三撃：a_m = 2sin(πm/k), 傾き Λ'(0⁺)=1 の証明
  branching-gate-error.md       第三撃：理想化の誤差評価（T→0 厳密, 有限T は O(T/k²)）
  nonsplittable-gate.md         Case III：capacity が係数依存の PDE-algebraic invariant
  nonsplittable-Ja.md           Case III：J(0)=1, J'(0)=-2π²/Γ(1/4)⁴（閉形式調査）
  nonsplittable-classification.md  Case III：一般 germ の J 分類（分離骨格＋Gamma結合則）
  branching-nonsplittable.md    多価×非分離：admittance が係数依存（調和点でのみ 2sin）
  open-questions.md / lean-design*.md
julia/       第一撃の横断積分（Julia）
  transverse_integral.jl            1 次元 Bessel 表示
  transverse_integral_direct_2d.jl  独立 2 次元検証
python/      数値実験（Python）
  split_gate_classification.py   第二撃：germ ライブラリ
  langevin_exit_time.py          第一撃：Langevin 平均脱出時間
  plot_exit_time.py              exit-time 検証図
  branching_gate_committor.py    第三撃：committor/capacity solver（Morse で検証）
  branching_gate_refine.py       第三撃：指数 p の精密化
  branching_valley_number.py     第三撃：谷数 k 依存
  check_environment.py
lean/        形式化のスケッチ（AbstractBridge / AsymptoticQuotient）
literature/  prior-art.md / claims-map.md
research-log/ 2026-07.md（実験ログ・数値結果）
data/ figures/  生成物（.csv/.png は .gitignore、コードで再生成）
```

## 再現手順

**Python 環境**

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python python/check_environment.py
```

**Julia 環境**（第一撃の横断積分に使用；Julia 1.12.x 推奨）

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. julia/check_environment.jl
```

**主な実験**

```bash
# 第一撃
julia --project=. julia/transverse_integral_direct_2d.jl   # 横断積分の独立2D検証
python python/langevin_exit_time.py --temps 0.08 0.07 0.06 0.055 0.05 --n-traj 600
python python/plot_exit_time.py

# 第二撃
python python/split_gate_classification.py

# 第三撃
python python/branching_gate_committor.py   # solver 検証（Morse）＋ monkey saddle
python python/branching_gate_refine.py       # 指数 p → 1 の確認
python python/branching_valley_number.py     # 谷数 k 依存
```

生成物は `data/*.csv`, `figures/*.png`（Git 管理外、上記コードで再生成）。

## 証拠の区分

- **厳密／既知理論一致**：横断積分の Bessel 表示と漸近；split-gate 公式の
  1D Kramers・多次元 Eyring–Kramers 再現（機械精度）；split gate の分類式；
  monkey saddle の `H=0`・`h_C=1/2`・Hessian=0 による split 不適用。
- **数値的に確定**：I(T) の独立 2D 検証；6 germ の `(λ⊥,m⊥)→(κ,r)`；Langevin の
  `(log(1/T))^{-1}` 署名；branching gate の `p=1`・谷数 k 依存。
- **予想／未決**：Newton-adapted capacity lower bound の厳密証明（theory.md §9）；
  branching-gate 定数 `J_k` の解析評価と谷数一般化；non-splittable gate（Case III）；
  非可逆 SGD への外挿。

## 主要な未解決問題

1. **厳密下界**：退化 gate の sharp capacity lower bound（Newton-adapted sublevel section）。
2. **branching 不変量**：committor germ / gate topology を含む capacity 分類不変量。
3. **非平衡拡張**：状態依存・異方的 noise（SGD 的）への外挿（theory.md Layer III）。

---

詳細な導出・数値結果・反証条件は各 `docs/` と `research-log/2026-07.md` を参照。
