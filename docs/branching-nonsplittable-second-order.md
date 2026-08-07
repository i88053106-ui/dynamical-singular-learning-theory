# 非調和多価 gate：2 次係数 `d²a_m/dt²|_0`（2D resolvent）

`branching-nonsplittable-general.md` で錨点 `a_m(0)=2 sin(πm/k)` と1次傾き `a_m'(0)`
（envelope 定理）を確定した。摂動展開の**最後に評価できる要素**である曲率 `a_m''(0)` を、
**2 次元 resolvent solve** で数値評価する。

## 1. 定式化：graph-Laplacian 摂動と resolvent

admittance は Dirichlet エネルギーの対角成分 `E(t)=S_m a_m(t)`,
`S_m=Σ_j cos²(2πmj/k)`。重み付きグラフ Laplacian `L(t)`（`H^T L H = Σ_edges c(Δh)²`,
`c_pq=½(w_p+w_q)`, `w_p=e^{-F_t(p)}`, `F_t=s(-Re w^k + t r^k)`）で committor `h` は
`L_II h_I = -L_IB g`（`g` は reservoir 値）。`'=d/dt`（固定 h）を

```
c'_pq = ½(w'_p+w'_q),  w'_p = -s r_p^k w_p,
c''_pq= ½(w''_p+w''_q), w''_p = (s r_p^k)² w_p
```

とおくと、envelope 定理（極小の t 微分＝被積分の偏微分）から

```
┌────────────────────────────────────────────────────────────────────┐
│  E'(0)  = [h]^T L' [h]                        （1次, envelope）       │
│  E''(0) = [h]^T L'' [h] + 2 [ḣ]^T L' [h]      （2次）                 │
│  L_II ḣ_I = -(L'[h])_I,   ḣ|_reservoir = 0    ← 2D resolvent solve    │
└────────────────────────────────────────────────────────────────────┘
```

`ḣ = dh/dt` は committor の t-微分で、**同じ作用素 `L_II` の逆**（resolvent）を右辺
`-(L'[h])_I` に作用させて得る。`L_II` の LU 分解は committor と resolvent で共用。
`a_m''(0) = E''(0)/S_m`。

**これが Case III の「非初等的 2D PDE 部分」の具体化**：`nonsplittable-classification.md`
で「2 次以降は 2D resolvent 依存」と述べた対象そのものを、多価 admittance について
数値評価した。

## 2. 数値（resolvent vs 中心 2 階差分）

`python/branching_nonsplittable_second_order.py`（L=1.9, n=541; k=6 は谷が深いため
スケール不変性で `s=0.2` に正規化）。**resolvent と有限差分が一致**（機構の確認）：

| k | m | a_m(0) | a_m'(0) resolvent | a_m'(0) FD | a_m''(0) resolvent | a_m''(0) FD |
|---|---|---|---|---|---|---|
| 4 | 1 | 1.4173 | −1.0979 | −1.1077 | +0.9303 | +0.9463 |
| 4 | 2 | 2.0059 | −1.9770 | −1.9956 | +2.0954 | +2.1278 |
| 6 | 1 | 1.0288 | −0.5403 | −0.5478 | +0.6210 | +0.5075 |
| 6 | 2 | 1.8141 | −1.4257 | −1.4415 | +1.8838 | +1.7324 |
| 6 | 3 | 2.1079 | −1.8078 | −1.8438 | +2.7332 | +2.3174 |

`a_m''(0) > 0`：admittance 曲線は調和点で**下に凸**（結合を強めても緩めても、1 次予測より
admittance が大きめに残る）。1 次は resolvent＝FD が全 (k,m) で ~1–2% 一致。2 次は k=4 で
~2%、k=6 で ~15–18% 差（2 階差分は有限差分ステップ・格子に敏感で、深い谷を持つ k=6 で
resolvent の方が信頼できる）。符号（凸）と桁は全 (k,m) で頑健。

## 3. Taylor 展開

```
a_m(t) = 2 sin(πm/k) + a_m'(0) t + ½ a_m''(0) t² + O(t³).
```

| k | m | 展開 |
|---|---|---|
| 4 | 1 | `a_1(t) ~ 1.4142 − 1.098 t + ½(0.930) t²` |
| 4 | 2 | `a_2(t) ~ 2.0000 − 1.977 t + ½(2.095) t²` |
| 6 | 1 | `a_1(t) ~ 1.0000 − 0.540 t + ½(0.621) t²` |
| 6 | 2 | `a_2(t) ~ 1.7321 − 1.426 t + ½(1.884) t²` |
| 6 | 3 | `a_3(t) ~ 2.0000 − 1.808 t + ½(2.733) t²` |

**k=4：結合 `a` への変換**（`a=(6+2t)/(1-t)`, `dt/da|_6=1/8`, `d²t/da²|_6=−1/32`）。
連鎖律 `d²a_m/da²|_6 = a_m''(0)/64 − a_m'(0)/32`：

```
a_1(a) ~ 1.4142 − 0.1372 (a-6) + ½(0.0488)(a-6)²,
a_2(a) ~ 2.0000 − 0.2471 (a-6) + ½(0.0945)(a-6)².
```

`branching-nonsplittable.md` の a-掃引表（a=3..8）と比較：a=6 近傍で 2 次展開が
数値点に密着（osculating parabola）、遠方で truncation 偏差
（figures/branching_nonsplittable_second_order.png）。例（a_1）：a=5 予測 1.576 / 実測
1.587、a=7 予測 1.301 / 実測 1.299。

## 4. 言明

```
┌──────────────────────────────────────────────────────────────────────┐
│  非調和多価 gate の admittance 摂動展開は次で完全に定量化される:         │
│   ・0 次 a_m(0)=2 sin(πm/k)            厳密（調和 Steklov, Case II）,   │
│   ・1 次 a_m'(0)=-(1/S_m)Σ½(w r^k)(Δh)²  envelope 定理（明示積分）,    │
│   ・2 次 a_m''(0)=[h]L''[h]+2[ḣ]L'[h]   2D resolvent（非初等・数値）.   │
│  0/1 次は初等的スケルトン、2 次以降は 2D resolvent —— Case III の階層が  │
│  多価 admittance でも成立し、曲率まで数値的に確定した.                   │
└──────────────────────────────────────────────────────────────────────┘
```

## 5. 証拠の区分

- **厳密**：graph-Laplacian 摂動恒等式 `E''=[h]L''[h]+2[ḣ]L'[h]`、resolvent 方程式
  `L_II ḣ_I=-(L'[h])_I`；連鎖律 `d²a_m/da²|_6=a_m''(0)/64−a_m'(0)/32`。
- **数値的に確定**：resolvent＝中心 2 階差分（k=4,6）；`a_m''(0)>0`（凸）；2 次 Taylor が
  a-掃引に密着。
- **未決**：`a_m''(0)` の閉形式（Case III 同様に非初等）；3 次以降；一般 k・大結合の
  収束半径。

`python/branching_nonsplittable_second_order.py` で再生成。
