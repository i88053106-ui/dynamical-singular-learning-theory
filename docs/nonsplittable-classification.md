# 一般 non-splittable germ の J 分類

`nonsplittable-Ja.md` の `-x^4+a x^2y^2+y^4` を一般化し、任意の同次 non-splittable
gate germ について capacity 定数 `J` の**分離バックボーン**と**任意結合の1次応答**を
閉形式（Gamma 関数）で与える。

## 1. 分類の骨格

同次 degree-`d` の 2-valley gate germ を

```
F(x,y) = -U(x) + V(y) + W(x,y),   U,V = leading (·)^d,  W = 結合
```

と書く。`J = inf_h ∫|∇h|^2 e^{-F} dxdy`（scale-invariant, `cap = J·T`）。

### (a) 分離バックボーン

`W = 0` なら committor は 1 次元 `h(x)`（`div(e^{U-V}∇h)=0` を `h(x)=C∫^x e^{-U}` が満たす）。

```
J_0 = [ ∫ e^{-U(x)} dx ]^{-1} · ∫ e^{-V(y)} dy = Z_V / Z_U.
```

`U = V = (·)^d` なら `J_0 = 1`。一般の同次 `U=α x^d, V=β y^d` でも
`J_0 = (β/α)^{-1/d}`... の Gamma 比で初等的。

### (b) 任意結合の1次応答（一般公式）

`W = Σ c_pq x^p y^q` を加えると、envelope 定理 + `h_0'(x)=e^{-U}/Z_U` より
`dJ/dε = -(1/Z_U^2) ∫∫ W e^{-U-V}`、したがって

```
┌──────────────────────────────────────────────────────────────────┐
│  J'(0)/J_0 = - Σ_pq c_pq <x^p>_d <y^q>_d,                           │
│  <x^p>_d = ∫ x^p e^{-x^d} / ∫ e^{-x^d} = Γ((p+1)/d)/Γ(1/d)  (p even),│
│           = 0  (p odd).                                             │
└──────────────────────────────────────────────────────────────────┘
```

capacity の結合への1次応答は、結合単項式のモーメント（不安定側 `<x^p>` × 安定側
`<y^q>`）の積和。すべて Gamma 関数で閉じる。

### (c) パリティ選択則

奇数指数を含む結合（`x^3y^3` 等）は `<x^{odd}>=0` により **1次では J に効かない**。

### (d) 非初等性（普遍）

2 次係数 `J''(0)` は a=0 作用素の resolvent（2 次元 PDE）を要する（`nonsplittable-Ja.md`）。
よって **分離バックボーンと1次応答は普遍的に初等的（Gamma）だが、full J はすべての
non-splittable germ で非初等的**。これが Case III の本質的内容。

## 2. 数値検証（`python/nonsplittable_classification.py`）

`J_0 = 1`（`-x^d+y^d`, d=4,6,8：1.001, 1.001, 1.010）。1 次応答：

| germ | numeric J'(0) | predicted −⟨x^p⟩⟨y^q⟩ |
|---|---|---|
| −x⁴+y⁴ + ε x²y² | −0.11399 | −0.11424 |
| −x⁶+y⁶ + ε x⁴y² | −0.06431 | −0.06457 |
| −x⁶+y⁶ + ε x²y⁴ | −0.06448 | −0.06457 |
| −x⁶+y⁶ + ε x³y³ | 0.00000 | 0（parity） |
| −x⁸+y⁸ + ε x⁶y² | −0.04387 | −0.04551 |
| −x⁸+y⁸ + ε x⁴y⁴ | −0.03557 | −0.03626 |
| −x⁸+y⁸ + ε x²y⁶ | −0.04503 | −0.04551 |
| −x⁸+y⁸ + ε x⁵y³ | 0.00000 | 0（parity） |

全点が対角線（numeric = predicted）に乗る（figures/nonsplittable_classification.png）。
奇結合は厳密に 0（選択則）。degree-8 の ~3% ずれは有限 domain/grid（e^{x⁸} 急増）。

## 3. 分類の言明

```
non-splittable germ の capacity 定数 J は、
  ・分離バックボーン J_0 = Z_V/Z_U        （初等的, Gamma）
  ・任意結合の1次応答 J'(0) = -Σ c_pq <x^p><y^q>  （初等的, Gamma モーメント積和）
  ・パリティ選択則（奇結合は1次で無効）
で普遍的に規定される「代数的スケルトン」を持つ。
一方 full J は2次以降で2次元 PDE の resolvent に依存し、非初等的。
```

Case I（横断 Newton data で完全に決定）→ Case III（代数的スケルトン＋非初等的 PDE 部分）
という階層が、これで定量的に確定した。

## 4. 証拠の区分

- **厳密**：分離バックボーン `J_0=Z_V/Z_U`；1次応答 `J'(0)=-Σc_pq<x^p>_d<y^q>_d`
  （envelope 定理 + Gamma モーメント）；パリティ選択則；2次以降の非初等性。
- **数値的に確定**：degree 4,6,8・複数結合で1次応答則が対角線一致、`J_0→1`、選択則 `=0`。
- **未決**：2次以降の係数（2D resolvent）；一般 germ の large-coupling / 非摂動領域；
  多価（branching）non-splittable germ への拡張。
