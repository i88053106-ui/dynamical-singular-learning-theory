# 動的特異学習理論

## Dynamical Singular Learning Theory

### A Singular Capacity Program for Learning Dynamics

---

## 一行要約

渡辺澄夫の特異学習理論が静的に分類する特異 Laplace 幾何を、確率的準安定力学の質量–capacity 構造へ接続し、とりわけ**特異遷移ゲートにおける small-noise capacity 漸近を分類する幾何学的不変量**を探究する。

本理論の基本図式は次である。

\[
\boxed{
\text{well geometry gives mass, gate geometry gives conductance, and time is their asymptotic ratio.}
}
\]

すなわち、

\[
\boxed{
\text{井戸の特異幾何は質量を与え、ゲートの特異幾何は伝導度を与える。学習時間は両者の漸近比として現れる。}
}
\]

可逆 Langevin 系では、準安定井戸の局所 Gibbs 質量

\[
Z_A(T)
\]

と transition capacity

\[
\operatorname{cap}_T(A,B)
\]

の比が exit-time scale を支配する。

象徴的には、

\[
Z_A(T)
\sim
C_A
e^{-L_A/T}
T^{\lambda_A}
\left(\log\frac1T\right)^{m_A-1},
\]

\[
\operatorname{cap}_T(A,B)
\sim
C_{A,B}
e^{-H_{A,B}/T}
T^{\kappa_{A,B}}
\left(\log\frac1T\right)^{r_{A,B}-1},
\]

なら、

\[
\mathbb E\tau_{A\to B}
\sim
C
T^{\lambda_A-\kappa_{A,B}}
\left(\log\frac1T\right)^{m_A-r_{A,B}}
\exp\left(
\frac{H_{A,B}-L_A}{T}
\right).
\]

ここで、

\[
\delta_{A\to B}
=
\lambda_A-\kappa_{A,B}
\]

という差そのものを新規性の中心とはしない。

これは mass–capacity asymptotics を整理する bookkeeping quantity である。

本理論の主要な数学問題は、

\[
\boxed{
\text{退化遷移ゲートの capacity asymptotics を何が分類するのか}
}
\]

である。

SLT における RLCT が singular Laplace integrals を分類するように、small-noise capacity を分類する特異不変量は存在するか。

それが potential germ の Newton data だけで決まるのか。

committor geometry を必要とするのか。

あるいは basin pair と gate network に本質的に依存するのか。

これが動的特異学習理論の数学的中心問題である。

---

# 0. 理論の三層構造と信頼性規約

本理論は三層からなる。

各層は異なる数学的問題を扱う。

---

## Layer I

### Singular capacity and metastable asymptotics

可逆 overdamped Langevin dynamics に対し、退化した transition gate の small-noise capacity asymptotics を研究する。

中心問題は、

\[
\operatorname{cap}_T(A,B)
\]

の前指数因子に現れる power exponent と logarithmic multiplicity を分類することである。

第一の定理候補は、

\[
\boxed{
\text{nondegenerate unstable direction}
+
\text{Newton-degenerate transverse gate}
}
\]

という split singular gate に限定する。

地位：定理化を直接狙う数学的核。

---

## Layer II

### Metastable state selection and Markov reduction

各準安定領域の局所 Gibbs 質量から局所自由エネルギーを定義し、長時間占有重みを記述する。

さらに低温 Langevin diffusion の metastable Markov reduction を考え、singular state weights と singular transition rates を持つ reduced network を構成する。

ここで本質的な問題は H-theorem 自体ではなく、

\[
\boxed{
\text{Langevin diffusion}
\longrightarrow
\text{metastable Markov network}
}
\]

という reduction の精度である。

地位：既知の entropy dissipation structure と singular mass/capacity asymptotics を接続する定理プログラム。

---

## Layer III

### Nonequilibrium singular metastability

SGD の状態依存・異方的 noise を含む非可逆小ノイズ力学へ拡張する。

一般には Gibbs measure も detailed balance も存在しない。

したがって RLCT を直接移植しない。

まず transition-time asymptotics を研究し、その sub-Arrhenius prefactor に再現可能な singular exponents が存在するかを問う。

地位：予想および研究課題。

---

# 1. 動機：静的特異幾何と学習時間の間の空白

特異学習理論では、階層的・過剰パラメータ化された統計モデルにおいて、最適分布を実現するパラメータ集合が一般に特異集合となる。

非負解析関数

\[
K(w)
\]

を平均 KL 情報量あるいは excess risk とすると、その零点集合

\[
W_0
=
\{w:K(w)=0\}
\]

は一般に滑らかな多様体ではない。

成分交差、高次退化、rank deficiency、parameter redundancy などが存在する。

SLT はこの局所特異幾何を real log canonical threshold

\[
\lambda
\]

と pole multiplicity

\[
m
\]

によって積分漸近へ変換する。

典型的には、

\[
\int e^{-K(w)/T}\varphi(w)\,dw
\sim
C
T^\lambda
\left(\log\frac1T\right)^{m-1}.
\]

この変換は、

\[
\boxed{
\text{singular geometry}
\longrightarrow
\text{mass scaling}
}
\]

である。

しかし学習には時間がある。

noisy gradient dynamics は、

- ある状態に長時間滞留し、
- plateau を形成し、
- 比較的短い時間に別の状態へ遷移し、
- 再び準定常状態を形成する、

という段階的挙動を示し得る。

問いたいのは、

> なぜその状態に長く滞留するのか。

> なぜある時点で別の状態へ移るのか。

> 特異幾何は transition clock にどのように現れるのか。

という問題である。

potential-theoretic metastability では、exit time は概念的に

\[
\frac{\text{well-side Gibbs weight}}
{\text{transition capacity}}
\]

という構造を持つ。

well-side integral は singular Laplace theory の対象である。

一方、capacity は equilibrium potential の Dirichlet energy を含む。

したがって本理論は、

\[
\boxed{
\text{SLT の Laplace geometry}
+
\text{metastability の capacity geometry}
}
\]

を一つの asymptotic calculus に置くことを目指す。

井戸側の singular Laplace contribution は、既存の singular Laplace asymptotics と接続しやすい。

これに対し、本理論が主要な数学的障害として扱うのは gate 側である。

---

# 2. 静的 SLT の最小構造

パラメータを

\[
w\in\mathbb R^d
\]

とし、非負解析関数 \(K(w)\) と滑らかな正の密度 \(\varphi(w)\,dw\) を考える。

局所ゼータ関数を

\[
\zeta(z)
=
\int K(w)^z\varphi(w)\,dw
\]

とする。

最大実部を持つ pole を

\[
-\lambda
\]

とし、その位数を \(m\) とする。

\(\lambda\) が RLCT である。

この pole data は sublevel volume

\[
V(\varepsilon)
=
\int_{K(w)<\varepsilon}
\varphi(w)\,dw
\]

の漸近

\[
V(\varepsilon)
\sim
C
\varepsilon^\lambda
|\log\varepsilon|^{m-1}
\]

および Laplace integral

\[
Z(T)
=
\int e^{-K(w)/T}\varphi(w)\,dw
\]

の漸近

\[
Z(T)
\sim
C'
T^\lambda
\left(\log\frac1T\right)^{m-1}
\]

を支配する。

SLT における本質は、

\[
\boxed{
\text{解析的特異 germ}
\longrightarrow
(\lambda,m)
\longrightarrow
\text{Laplace scaling}
}
\]

という分類である。

DSLT の中心的な問いは、

> capacity に対して、これに対応する分類理論を構成できるか。

である。

---

# 3. Langevin dynamics と capacity の規約

可逆 overdamped Langevin dynamics

\[
dW_t
=
-\nabla L(W_t)\,dt
+
\sqrt{2T}\,dB_t
\]

を考える。

\(L\) は confinement 条件を満たす解析的 potential とする。

非正規化 Gibbs measure を

\[
\mu_T(dw)
=
e^{-L(w)/T}\,dw
\]

とする。

Dirichlet form を

\[
\boxed{
\mathcal E_T(h,h)
=
T\int
|\nabla h(w)|^2
e^{-L(w)/T}\,dw
}
\]

と定義する。

互いに素な集合 \(A,B\) に対し、

\[
\operatorname{cap}_T(A,B)
=
\inf_h
\mathcal E_T(h,h),
\]

ただし

\[
h|_A=1,
\qquad
h|_B=0.
\]

minimizer \(h_{A,B}\) は equilibrium potential、すなわち committor に対応する。

適切な metastable regime と初期分布の下では、mean transition time は概念的に

\[
\boxed{
\mathbb E\tau_{A\to B}
\sim
\frac{\text{well-side Gibbs contribution}}
{\operatorname{cap}_T(A,B)}
}
\]

という mass–capacity structure を持つ。

重要なのは、分子と分母の数学的性質が非対称なことである。

well-side term は Laplace integral に近い。

capacity は

\[
T\int
|\nabla h_{A,B}|^2
e^{-L/T}\,dw
\]

という PDE variational quantityである。

---

# 4. 新規性の焦点：singular transition capacity

本理論では、

\[
\boxed{
\text{well-side degeneracy is not the primary novelty target}
}
\]

という立場をとる。

孤立した退化井戸では、well-side Gibbs integral を singular Laplace asymptotics によって評価することが主要な構造となる。

井戸側の RLCT pole data は依然として重要であるが、その役割は metastable mass の分類にある。

これに対し、本理論が数学的新規性の中心として扱うのは、

\[
\boxed{
\text{singular transition capacity}
}
\]

である。

また、

\[
\delta_{A\to B}
=
\lambda_A-\kappa_{A,B}
\]

を metastable singularity gap と呼ぶ。

ただしこれは新しい「差の原理」を主張する概念ではなく、well mass と capacity の asymptotic bookkeeping を行う edge quantity として用いる。

---

# 5. 井戸側の pole data

準安定領域 \(A\) に対して、

\[
L_A
=
\inf_{w\in A}L(w)
\]

とする。

局所 Gibbs mass を

\[
Z_A(T)
=
\int_Ae^{-L(w)/T}\,dw
\]

と定義する。

低温極限で、

\[
\boxed{
Z_A(T)
\sim
C_A
e^{-L_A/T}
T^{\lambda_A}
\left(\log\frac1T\right)^{m_A-1}
}
\]

とする。

\[
(\lambda_A,m_A)
\]

を well-side pole data と呼ぶ。

一つの解析的 germ が basin mass を支配する場合、これはその local zeta pole data と一致し得る。

複数の最低値成分が競合する場合には、basin integral 全体の leading asymptotics として operational に定義する。

本理論の立場は明確である。

\[
\boxed{
\text{well pole data classify mass, but do not by themselves classify clocks.}
}
\]

---

# 6. Capacity pole data と metastable singularity gap

transition capacity に対し、

\[
\operatorname{cap}_T(A,B)
\sim
C_{A,B}
e^{-H_{A,B}/T}
T^{\kappa_{A,B}}
\left(\log\frac1T\right)^{r_{A,B}-1}
\]

という sharp asymptotic が存在すると仮定する。

\[
(\kappa_{A,B},r_{A,B})
\]

を capacity asymptotic data と呼ぶ。

これは最初から saddle RLCT と同一視しない。

capacity は、

- local potential germ,
- committor geometry,
- stable/unstable splitting,
- gate topology,
- basin pair,
- competing channels,

に依存する可能性がある。

power exponent の差を

\[
\delta_{A\to B}
=
\lambda_A-\kappa_{A,B}
\]

とする。

logarithmic exponent の差を

\[
q_{A\to B}
=
m_A-r_{A,B}
\]

とする。

この組

\[
(\delta_{A\to B},q_{A\to B})
\]

を transition-edge asymptotic data とみなす。

これは状態不変量ではない。

direction-dependent な edge quantity である。

---

# 7. Singular Eyring–Kramers accounting principle

mass–capacity relation と二つの漸近を組み合わせれば、形式的に

\[
\boxed{
\mathbb E\tau_{A\to B}
\sim
C
T^{\delta_{A\to B}}
\left(\log\frac1T\right)^{q_{A\to B}}
\exp\left(
\frac{\Delta H_{A\to B}}{T}
\right)
}
\]

を得る。

ここで、

\[
\Delta H_{A\to B}
=
H_{A,B}-L_A.
\]

対数表示では、

\[
\log\mathbb E\tau_{A\to B}
=
\frac{\Delta H_{A\to B}}{T}
+
\delta_{A\to B}\log T
+
q_{A\to B}\log\log\frac1T
+
O(1).
\]

本理論ではこれを新定理とは呼ばない。

これは accounting principle である。

本当の問題は、

\[
\boxed{
(\kappa,r)\text{ をどの幾何が決めるのか}
}
\]

である。

---

# 8. 第一の定理候補：Newton-degenerate transverse gate

第一標的として、三次元解析的 potential

\[
\boxed{
F(x,y,z)
=
\frac{(x^2-1)^2}{4}
+
x^2(y^2+z^2)
+
y^6+y^2z^2+z^6
}
\]

を考える。

横断 singular germ を

\[
K(y,z)
=
y^6+y^2z^2+z^6
\]

とする。

critical points は、

\[
(-1,0,0),
\qquad
(0,0,0),
\qquad
(1,0,0)
\]

のみである。

\((\pm1,0,0)\) は Morse minima。

原点では Hessian eigenvalues は

\[
(-1,0,0).
\]

left well から right well への任意の path は \(x=0\) を横切る。

その断面では、

\[
F(0,y,z)
=
\frac14+K(y,z)
\geq
\frac14.
\]

equality は \((y,z)=(0,0)\) のみである。

したがって communication height は

\[
H=\frac14
\]

であり、原点が unique gate である。

原点近傍では、

\[
F-\frac14
=
-\frac{x^2}{2}
+
x^2(y^2+z^2)
+
\frac{x^4}{4}
+
K(y,z).
\]

さらに局所解析座標

\[
X
=
x
\sqrt{
1-2(y^2+z^2)-\frac{x^2}{2}
}
\]

を用いれば、原点近傍で

\[
\boxed{
F-\frac14
=
-\frac{X^2}{2}
+
K(y,z)
}
\]

と exact split normal form にできる。

したがってこの gate は、

\[
\boxed{
\text{nondegenerate unstable direction}
+
\text{Newton-degenerate transverse geometry}
}
\]

を持つ。

横断 Laplace integral

\[
I(T)
=
\int_{\mathbb R^2}
e^{-K(y,z)/T}\,dy\,dz
\]

は、

\[
\boxed{
I(T)
\sim
C_K
T^{1/2}
\log\frac1T
}
\]

という pole multiplicity two 型漸近を持つとする。

候補定数は、

\[
C_K=\frac{\sqrt\pi}{3}
\]

である。

これは別途厳密に証明する。

capacity の予測は、

\[
\boxed{
\operatorname{cap}_T(A,B)
\sim
\frac{1}{3\sqrt2}
T\log\frac1T
\exp\left(-\frac{1}{4T}\right)
}
\]

である。

well minima は Morse で Hessian determinant は \(8\)。

したがって exit-time prediction は、

\[
\boxed{
\mathbb E\tau_{-\to+}
\sim
3\sqrt2\,\pi^{3/2}
\frac{T^{1/2}}{\log(1/T)}
\exp\left(\frac{1}{4T}\right)
}
\]

である。

この

\[
\boxed{
(\log(1/T))^{-1}
}
\]

という負の logarithmic prefactor exponent が、第一標的の観測署名である。

---

# 9. Newton-adapted capacity lower bound

第一定理候補の中心的 proof obligation は、既存の fiberwise capacity lower-bound architecture を Newton geometry に適応させることである。

横断 Euclidean ball の代わりに、

\[
\boxed{
D_T
=
\left\{
(y,z):
K(y,z)
<
M T\log\frac1T
\right\}
}
\]

という Newton-adapted sublevel section を用いる。

Dirichlet energy を局所領域へ制限し、

\[
|\nabla h|^2
\geq
|\partial_X h|^2
\]

とする。

固定した transverse coordinate \((y,z)\) ごとに unstable direction \(X\) 上の一変数変分問題へ還元する。

この fiberwise reduction では transverse domain の側面形状は主項に直接現れない。

必要なのは、

1. 両端の lids 上で equilibrium potential が \(0\) または \(1\) に一様に近いこと、
2. \(D_T\) が transverse Gibbs mass の \(1-o(1)\) を捕捉すること、
3. 座標変換の Jacobian と Dirichlet metric correction が \(1+o(1)\) であること、

である。

mass capture は、

\[
I(T)
\sim
C_KT^{1/2}\log(1/T)
\]

を用いれば短く示せる。

\(D_T^c\) 上では、

\[
e^{-K/T}
\leq
T^{M/2}e^{-K/(2T)}.
\]

したがって、

\[
\int_{D_T^c}e^{-K/T}
\leq
T^{M/2}I(2T)
=
o(I(T)).
\]

よって、

\[
\boxed{
\int_{D_T}e^{-K/T}
=
(1-o(1))I(T).
}
\]

第一論文の focused extension は、

\[
\boxed{
\text{Euclidean transverse sections}
\longrightarrow
\text{Newton-adapted singular sublevel sections}
}
\]

への sharp capacity lower-bound extension である。

---

# 10. Singular capacity の三分類

理論構造と第一モデルの解析を踏まえ、singular gates を三つに分ける。

---

## Case I: split singular gate

局所座標で、

\[
F-H
=
-u(s)
+
K(v)
+
\text{higher-order terms}
\]

と stable/unstable structure が分離できる。

特に unstable direction が非退化なら、

\[
u(s)=\frac{\mu s^2}{2}.
\]

この場合、capacity は transverse singular Laplace integral に還元される可能性が高い。

第一モデルはこの class に属する。

ここでは Newton polyhedron、resolution of singularities、RLCT-type pole data が capacity prefactor を分類する候補となる。

---

## Case II: branching gate

gate 近傍で三つ以上の downhill components が存在する。

committor の局所値が単純な \(0/1\) boundary layer では決まらない。

複数の valley 間の harmonic splitting が必要になる。

この場合、

\[
\text{potential germ alone}
\]

では capacity data を決められない可能性がある。

gate network と local committor geometry が必要になる。

---

## Case III: non-splittable singular gate

unstable direction と singular kernel geometry が局所的に分離できない。

単一の nondegenerate unstable quadratic mode を split off できない。

この場合、

\[
|\nabla h|^2e^{-F/T}
\]

の singular perturbation structure を potential germ と同時に解析する必要がある。

ここで初めて、

\[
\boxed{
\text{PDE–algebraic capacity invariant}
}
\]

が必要になる可能性がある。

本理論の最終数学目標は Case II–III の分類である。

---

# 10′. 三分類の到達点：確定した capacity・admittance 構造

§10 の三分類に対し、明示モデルの解析と数値検証で確定した構造を総合する。
到達点は一つの階層に収まる。

\[
\boxed{
\begin{aligned}
&\textbf{Case I : capacity は横断 Newton pole data で決まる（初等的）.}\\
&\textbf{Case II : branching network の admittance は調和性ゆえ普遍 } a_m=2\sin(\pi m/k)\ \textbf{（初等的）.}\\
&\textbf{Case III: network は残るが admittance・capacity は係数依存の PDE-algebraic 量.}\\
&\qquad\textbf{普遍性は調和点でのみ. 代数骨格は初等的, full 値は 2D resolvent（非初等）.}
\end{aligned}
}
\]

各 case の committor 構造と、capacity/admittance を決める幾何、その初等性：

| case | committor | capacity/admittance を決めるもの | 初等性 |
|---|---|---|---|
| I split | fiberwise 1 次元 | 横断 Newton pole data \((\lambda_\perp,m_\perp)\) | 初等（Newton/RLCT）|
| II branching | 2 次元 harmonic | Steklov 固有値 \(|m|\) → \(a_m=2\sin(\pi m/k)\) | 初等（普遍・係数非依存）|
| III non-splittable | 2 次元 非分離 | 分離骨格（Gamma）＋ 2D resolvent | 骨格 初等 / full 非初等 |
| 多価 × 非分離 | 2 次元 harmonic 網＋非調和 | 調和点で \(2\sin\)、周囲は resolvent 摂動 | 0/1 次 初等 / 2 次〜 非初等 |

---

## Case I（split gate）：横断 Newton pole data による分類 — 確定

第一撃 model \(F=(x^2-1)^2/4+x^2(y^2+z^2)+y^6+y^2z^2+z^6\)、gate \(H=1/4\)、横断 germ
\(K=y^6+y^2z^2+z^6\)。exact split normal form \(F-1/4=-X^2/2+K(y,z)\)。

- 横断積分 \(I(T)=(\sqrt\pi/3)\sqrt T\,[\log(1/T)+\gamma+6\log2]\)（厳密 Bessel 表示
  \(I=(4/3)\sqrt T\int e^{-x^2}K_0(2\sqrt T x^3)dx\) と独立 2D が相対差 \(10^{-10}\) で一致）。
- capacity \(\operatorname{cap}_T\sim(1/3\sqrt2)\,T[\log(1/T)+\gamma+6\log2]\,e^{-1/4T}\),
  exit \(\mathbb E\tau\sim 3\sqrt2\,\pi^{3/2}\,\sqrt T/\log(1/T)\,e^{1/4T}\)（Eyring–Kramers を機械精度で再現）。
- Langevin 7 温度で負 log prefactor 署名を確認（\(C_{\rm eff}=2.055\pm0.025\)）。

第二撃（一般 split gate）：横断 germ の Laplace pole data から

\[
\boxed{(\kappa,r)=(\lambda_\perp+\tfrac12,\ m_\perp)},\qquad
\lambda_\perp=1/d\ (\text{Newton 距離}),\quad m_\perp\in\{1,2\}\ (\text{辺 or 頂点}).
\]

6 germ で検証。決定的対比：\(y^4+y^2z^2+z^4\)（中央項が辺上 → \(m=1\), log なし）vs
\(y^6+y^2z^2+z^6\)（頂点 → \(m=2\), log あり）。log multiplicity は Hessian でなく
Newton 頂点条件が支配（**予想 B 確定**）。文書：`capacity-derivation.md`,
`split-gate-classification.md`。

---

## Case II（branching gate）：普遍 Steklov admittance — 確定

monkey saddle \(F=r^6/6-(x^3-3xy^2)\) および \(\mathrm{Re}(w^k)\) 族。split 分解は原理的に
不適用（不安定曲率 \(\mu=0\)）。

- branching 署名 \(h_C=1/2\)（committor が 2D harmonic の直接証拠）、capacity 指数 \(p=1\)
  （log なし、degree-3 homogeneous scaling）、谷数 k で power は不変・定数は単調減少。
- circulant k 端子ネットワーク：\(R_{0d}=(1/k)\sum_m 4\sin^2(\pi md/k)/a_m\)。
- **普遍 admittance**（単位円 Steklov 固有値 \(|m|\)、conformal \(W=w^k\) と分数次数ベッセル
  \(I_{|m|/k}\) 由来）：

\[
\boxed{a_m=2\sin(\pi m/k)},\qquad \Lambda(q)=2\sin(q/2),\ \Lambda'(0^+)=1,\qquad
J_k=\tfrac{k}{2}\tan\tfrac{\pi}{2k}\ \xrightarrow{k\to\infty}\ \tfrac\pi4.
\]

- 有限 T 補正 \(a_m(k,T)=2\sin(\pi m/k)(1-c_m(k)T+O(T^2))\)、\(c_m\sim c_m/k^2\)、\(T\to0\) で厳密。

これらは germ が**調和**（\(\mathrm{Re}\,w^k\)）であることに由来し、**係数非依存**。文書：
`branching-gate.md`, `-Jk.md`, `-bs.md`, `-largek.md`, `-slope.md`, `-error.md`。

---

## Case III（non-splittable gate）：係数依存 PDE-algebraic capacity — 確定

同次 2-valley germ \(F=-U(x)+V(y)+W(x,y)\)、\(J=\inf_h\int|\nabla h|^2 e^{-F}\)
（scale-invariant, \(\operatorname{cap}=J\cdot T\)）。

- 分離バックボーン（\(W=0\)、committor 1 次元）：\(J_0=Z_V/Z_U=\int e^{-V}/\int e^{-U}\)、
  \(U=V=(\cdot)^d\) で \(J_0=1\)。
- 任意結合の1次応答（envelope 定理 + Gamma モーメント）：

\[
\boxed{J'(0)/J_0=-\sum_{pq}c_{pq}\langle x^p\rangle_d\langle y^q\rangle_d},\qquad
\langle x^p\rangle_d=\frac{\Gamma((p+1)/d)}{\Gamma(1/d)}\ (p\ \text{偶}),\ 0\ (p\ \text{奇}).
\]

パリティ選択則（奇結合は1次で無効）。degree 4,6,8・複数結合で検証。pure germ
\(-x^4+a x^2y^2+y^4\)：\(J(0)=1\), \(J'(0)=-2\pi^2/\Gamma(1/4)^4\approx-0.1142\)。

- **full \(J\) は2次以降 2D resolvent 依存（非初等）**。よって Case III の \(J\) は
  「代数的スケルトン（初等的 Gamma）＋非初等的 PDE 部分」（**予想 C の non-splittable 部分を
  定量化**）。Newton polygon \(\{(4,0),(2,2),(0,4)\}\) 一定でも \(J(a)\) は係数で変わる
  （Case I との決定的対比）。文書：`nonsplittable-gate.md`, `-Ja.md`, `-classification.md`。

---

## 多価 × 非分離：Case II と Case III の合流 — 確定

族 \(P_a=-x^4+ax^2y^2-y^4\)（\(a>2\) で 4 谷、\(\Delta P_a=(2a-12)r^2\) より調和は \(a=6\) の
\(-\mathrm{Re}\,w^4\) のみ）。再パラメータ \(P_{k,t}=-\mathrm{Re}(w^k)+t\,(x^2+y^2)^{k/2}\)
（k 偶、\(t=0\) 調和）。等方項 \(t\,r^k\) は谷を動かさないので circulant k 端子ネットは全 t で不変。

- **network は保たれ、admittance は係数依存**：普遍 \(2\sin(\pi m/k)\) は調和点 \(t=0\) のみ。
  スケール不変性で \(a_m\) は**単一変数** \(t\) の関数。k=4 で \(a=(6+2t)/(1-t)\)。
- **調和点まわりの摂動展開**（\(S_m=\sum_j\cos^2(2\pi mj/k)\)、\(h\)=調和 committor、\(L(t)\)=重み
  付きグラフ Laplacian）：

\[
\boxed{
\begin{aligned}
&\text{0 次}: a_m(0)=2\sin(\pi m/k) &&(\text{厳密, Steklov}),\\
&\text{1 次}: a_m'(0)=-\tfrac1{S_m}\!\sum_{\rm edges}\!\tfrac12(w_p r_p^k+w_q r_q^k)(\Delta h)^2 &&(\text{envelope, 明示積分}),\\
&\text{2 次}: a_m''(0)=[h]^\top L''[h]+2[\dot h]^\top L'[h],\ \ L_{II}\dot h_I=-(L'[h])_I &&(\text{2D resolvent, 非初等}).
\end{aligned}
}
\]

\(a_m''>0\)（凸）。全て resolvent＝中心差分で検証（k=4,6）。k=4 の連鎖律：
\(da_m/da|_6=\tfrac18 a_m'(0)\)（\(=-0.137,-0.247\)）、\(d^2a_m/da^2|_6=a_m''(0)/64-a_m'(0)/32\)。

- **収束半径**：\(a_m(t)\) は \(t\in(-1,1)\)（\(a\in(2,\infty)\)）で解析的、両端で degree-4
  同次のまま退化——

\[
t\to-1:\ P\to-(x^2-y^2)^2\ (\text{井戸合体},\ a_m\sim(a-2)^{-1/2});\quad
t\to+1:\ P\to 8x^2y^2\ (\text{Newton 退化},\ a_m\sim a^{-\gamma_+},\ \gamma_+\approx0.15).
\]

いずれも非整数冪の分岐点なので \(\boxed{R=1}\)（実用は \(|t|\lesssim0.3\)）。文書：
`branching-nonsplittable.md`, `-general.md`, `-second-order.md`, `-radius.md`。

\[
\boxed{
\begin{aligned}
&\textbf{多価 non-splittable gate の大域像：}\\
&\text{・network は circulant（Case II 構造）を全域 }|t|<1\text{ で保つ,}\\
&\text{・admittance }a_m(t)\text{ は調和点 }t=0\text{ で普遍 }2\sin(\pi m/k)\text{、周囲で係数依存の解析関数（Case III）,}\\
&\text{・展開は 0 次(Steklov 厳密)・1 次(envelope)・2 次(2D resolvent) で曲率まで確定, 収束半径 }R=1,\\
&\text{・両端 }t=\pm1\text{ は Newton 退化（capacity の分岐点）.}\\
&\textbf{Case II（普遍）と Case III（係数依存）は「network 共通・普遍性は調和性由来」で統一される.}
\end{aligned}
}
\]

---

# 11. 最小次元について

separable analytic non-branching one-gate class に限定する。

transverse dimension が一である場合、非負一変数 analytic germ は leading order で

\[
cy^{2k}
\]

型となり、その Laplace integral は pure power scaling を持つ。

したがって pole-multiplicity 型の leading logarithmic transverse factorを生むには transverse dimension が少なくとも二必要である。

よって、

\[
\boxed{
\text{separable analytic non-branching class では ambient dimension }d\geq3
}
\]

が必要である。

これは一般の二次元 singular capacity に log factor が存在しないことを意味しない。

branching、non-splittable、nonanalytic geometry は別問題である。

---

# 12. 準安定 basin の局所自由エネルギー

各 basin \(A_i\) に対し、

\[
Z_i(T)
=
\int_{A_i}e^{-L/T}\,dw
\]

とする。

局所自由エネルギーを

\[
\mathcal F_i^{\rm loc}(T)
=
-T\log Z_i(T)
\]

と定義する。

もし、

\[
Z_i(T)
\sim
C_i
e^{-L_i/T}
T^{\lambda_i}
\left(\log\frac1T\right)^{m_i-1},
\]

なら、

\[
\boxed{
\mathcal F_i^{\rm loc}(T)
=
L_i
+
\lambda_iT\log\frac1T
-
(m_i-1)T\log\log\frac1T
-
T\log C_i
+
o(T)
}
\]

である。

leading singular correction は、

\[
L_i
+
\lambda_iT\log\frac1T.
\]

ただし multiplicity が競合する場合には、

\[
-(m_i-1)T\log\log\frac1T
\]

も状態選択に寄与する。

この量は新しい独立熱力学量ではない。

局所 Gibbs mass の漸近を書き直したものである。

---

# 13. 占有重みと transition clock は異なる

\[
Z_i(T)
\]

と

\[
\mathcal F_i^{\rm loc}(T)
\]

は basin の equilibrium occupation weight を記述する。

一方、

\[
\operatorname{cap}_T(A_i,A_j)
\]

と

\[
\mathbb E\tau_{i\to j}
\]

は transition graph と transition clock を記述する。

したがって、

\[
\mathcal F_j^{\rm loc}
<
\mathcal F_i^{\rm loc}
\]

であっても、\(i\to j\) が直接または即座に起こるとは限らない。

基本図式は、

\[
\boxed{
\text{occupation weights}
\leftrightarrow
\text{local free energies}
}
\]

\[
\boxed{
\text{transition routes and clocks}
\leftrightarrow
\text{capacity network}
}
\]

である。

---

# 14. Coarse-grained free energy

metastable basin の占有確率を

\[
p_i(t)
=
\Pr(W_t\in A_i)
\]

とする。

time-scale separation の下で reduced Markov chain

\[
\dot p=pQ
\]

が得られる regime を考える。

reduced equilibrium weight を

\[
\pi_i
=
\frac{Z_i(T)}
{Z_{\rm cg}},
\qquad
Z_{\rm cg}
=
\sum_jZ_j(T)
\]

とする。

coarse-grained free energyを、

\[
\mathcal G_{\rm cg}(p,T)
=
\sum_ip_i\mathcal F_i^{\rm loc}(T)
+
T\sum_ip_i\log p_i
\]

と定義する。

すると、

\[
\boxed{
\mathcal G_{\rm cg}(p,T)
=
T D_{\rm KL}(p\Vert\pi)
-
T\log Z_{\rm cg}.
}
\]

したがって reduced reversible chain 上では relative entropy dissipation により、

\[
\frac{d}{dt}
\mathcal G_{\rm cg}(p(t),T)
\leq0.
\]

新しい問題は H-theorem 自体ではない。

問題は singular state weights と singular transition rates を持つ Markov reduction を、連続 diffusion からどの精度で導出できるかである。

---

# 15. Metastable Markov network

各 edge \(i\to j\) に、

\[
\operatorname{cap}_T(A_i,A_j)
\sim
C_{ij}
e^{-H_{ij}/T}
T^{\kappa_{ij}}
\left(\log\frac1T\right)^{r_{ij}-1}
\]

を対応させる。

各 node \(i\) には、

\[
Z_i(T)
\sim
C_i
e^{-L_i/T}
T^{\lambda_i}
\left(\log\frac1T\right)^{m_i-1}
\]

を対応させる。

reduced rate は象徴的に、

\[
q_{ij}(T)
\asymp
\frac{\operatorname{cap}_T(A_i,A_j)}
{Z_i(T)}.
\]

したがって transition network は、

\[
\boxed{
\text{node mass data}
+
\text{edge capacity data}
}
\]

によって記述される。

DSLT では、learning development を singular state sequence と同一視しない。

まず metastable network が存在するかを確認する。

その後、algebraic strata と network states の対応を研究する。

---

# 16. SGD の小ノイズ力学

SGD の連続時間近似を象徴的に、

\[
dW_t
=
b(W_t)\,dt
+
\sqrt{\varepsilon}\sigma(W_t)\,dB_t
\]

と書く。

\[
\Sigma(w)
=
\sigma(w)\sigma(w)^\top
\]

は一般に状態依存・異方的である。

drift と diffusion が Gibbs reversible structure を満たす保証はない。

したがって、

\[
e^{-L/T}
\]

を stationary measure と仮定しない。

また、

\[
L+\lambda T\log(1/T)
\]

をSGDへ直接移植しない。

---

# 17. 非可逆系では transition-time asymptotics から始める

まず、

\[
\mathbb E\tau_{A\to B}
=
e^{\Delta V_{A\to B}/\varepsilon}
a_{A,B}(\varepsilon)
\]

と書く。

\[
\Delta V_{A\to B}
\]

は quasipotential barrier に対応する。

次に、

\[
\gamma_{A\to B}
=
\lim_{\varepsilon\to0}
\frac{
\log a_{A,B}(\varepsilon)
}{
\log\varepsilon
}
\]

が存在する条件を問う。

さらに polyhomogeneous expansion

\[
a_{A,B}(\varepsilon)
\sim
C
\varepsilon^\gamma
\left(\log\frac1\varepsilon\right)^q
\]

が成立するかを研究する。

well–gate decomposition

\[
\gamma
=
\chi_A-\chi_{A,B}^{\rm gate}
\]

は最初から仮定しない。

まず transition prefactor が regular variation structure を持つかを問う。

---

# 18. Local learning coefficient との関係

empirical local learning coefficientを、

\[
\hat\lambda(t)
\]

とする。

少なくとも三つのレベルを区別する。

1. analytic germ に付随する asymptotic pole data、
2. finite-scale effective local exponent、
3. numerical estimator \(\hat\lambda(t)\)。

\[
\hat\lambda(t)
\]

を true RLCT と直接同一視しない。

また、

\[
\hat\lambda(t)
\]

が transition-edge data

\[
(\delta,q)
\]

を直接与えるとも仮定しない。

LLC の役割は、

> trajectory が異なる singular neighborhoods を横断している可能性を示す local geometric probe

である。

DSLT の強い主張候補は、

\[
\boxed{
\text{local singular geometry}
\longrightarrow
\text{transition prefactor data}
}
\]

という clock law である。

---

# 19. 中心予想

## 予想 A

### Logarithmic capacity at a Newton-degenerate saddle

第一モデルについて、

\[
\operatorname{cap}_T(A,B)
\sim
C
T\log\frac1T
e^{-1/(4T)}
\]

が成立する。

従って、

\[
\mathbb E\tau
\sim
C'
\frac{T^{1/2}}{\log(1/T)}
e^{1/(4T)}.
\]

**状態：第一撃で確定**（§10′ Case I）。\(\operatorname{cap}_T\sim(1/3\sqrt2)T[\log(1/T)+\gamma+6\log2]e^{-1/4T}\)、
Langevin 7 温度で負 log prefactor 署名 \(C_{\rm eff}=2.055\pm0.025\)。

---

## 予想 B

### Split singular gates are classified by transverse pole data

局所的に、

\[
F-H
=
-\frac{\mu s^2}{2}
+
K(v)
\]

へ還元される one-gate model について、\(K\) の singular Laplace pole data

\[
(\lambda_\perp,m_\perp)
\]

が capacity prefactor の power/log data を決める。

**状態：第二撃で確定**（§10′ Case I）。\((\kappa,r)=(\lambda_\perp+\tfrac12,m_\perp)\)、
\(\lambda_\perp=1/d\)、\(m_\perp\in\{1,2\}\)（辺 or 頂点）、6 germ で検証。

---

## 予想 C

### Non-splittable and branching gates require additional geometry

unstable mode を singular transverse geometry から分離できない場合、または branching gate の場合、

\[
\text{potential pole data alone}
\]

ではcapacity asymptoticsを分類できない。

committor germ、local harmonic splitting、gate topology の追加情報が必要になる。

**状態：Case II/III で定量化**（§10′）。branching（調和）は普遍 admittance \(a_m=2\sin(\pi m/k)\)
（Steklov \(|m|\)、係数非依存）；non-splittable は係数依存の PDE-algebraic 量で、分離骨格
（Gamma）＋ 2D resolvent。多価×非分離では両者が「network 共通・普遍性は調和点のみ」で合流し、
調和点まわりの摂動展開（0/1/2 次）と収束半径 \(R=1\) まで確定。潜在情報（committor geometry・
gate topology）が本質的に必要という予想 C の主張は、Case III で具体化された。

---

## 予想 D

### Occupation weights and transition clocks are governed by different structures

basin occupation weights は local Gibbs masses が支配する。

transition clocks は capacity network が支配する。

representation development の順序は local free-energy rankingだけでは決まらない。

---

## 予想 E

### Nonequilibrium transition exponents may exist for SGD-like diffusion

適切な small-noise regime では、

\[
\log\mathbb E\tau
=
\frac{\Delta V}{\varepsilon}
+
\gamma\log\varepsilon
+
q\log\log\frac1\varepsilon
+
O(1)
\]

という expansion が存在するクラスがある。

その分類不変量は可逆極限でwell/capacity pole dataへ還元される可能性がある。

---

# 20. Grokking 仮説

grokking plateau を metastable residence としてモデル化できる regime を考える。

memorization-like basinを \(M\)、generalization-like basinを \(G\) とする。

DSLTは、

> grokking は basin transition である

という一般論を独自主張とはしない。

より強い問いは、

\[
\boxed{
\text{noise scale に対する transition clock の sub-Arrhenius slope を予測できるか}
}
\]

である。

象徴的には、

\[
\log\tau_{\rm grok}
=
\frac{\Delta V}{\varepsilon}
+
\gamma\log\varepsilon
+
q\log\log\frac1\varepsilon
+
O(1).
\]

learning rate、batch size、explicit noise、weight decay は、

- quasipotential barrier、
- local metastable mass、
- transition gate geometry、

を別々に変形する可能性がある。

ただし実際のSGDにmetastable first-exit pictureが成立するかは独立に検証する。

---

# 21. Saddle-to-saddle learning との接続

deep linear network、matrix factorization、low-rank learningでは、段階的なrank/feature acquisitionが観測され得る。

本理論はこれを直ちにRLCTの単調変化とは解釈しない。

必要なのは、

\[
S_r
\]

という algebraic rank strata と、

\[
A_i
\]

という dynamical metastable basins の対応を調べることである。

trajectory

\[
A_{i_1}
\to
A_{i_2}
\to
\cdots
\]

に対し、各node mass dataとedge capacity dataを推定する。

対応が成立するなら、

\[
\boxed{
\text{feature acquisition}
\longleftrightarrow
\text{motion on a singular metastable network}
}
\]

という力学像が得られる。

---

# 22. Flat minima 理論との違い

flat minima theory はHessian eigenvaluesやcurvatureを中心に議論されることが多い。

しかし特異点ではHessianは局所積分体積を完全には特徴づけない。

例えば、

\[
K(y,z)
=
y^6+y^2z^2+z^6
\]

ではHessianは零行列である。

しかしcapacityに必要な transverse Laplace mass は、

\[
T^{1/2}\log(1/T)
\]

という特定の scaling を持つ。

zero eigenvalue count だけからこのlogarithmic multiplicityは復元できない。

したがって本理論の主張は、

> flatness がescape timeを決める

ではない。

より正確には、

\[
\boxed{
\text{singular mass geometry and singular conductance geometry jointly determine metastable clocks}
}
\]

である。

---

# 23. 最小検証モデル

## A. Singular Laplace toy models

\[
x^2y^2
\]

や、

\[
x^6+x^2y^2+y^6
\]

を用いる。

目的はmetastabilityではない。

\[
Z(T)
\sim
CT^\lambda
(\log(1/T))^{m-1}
\]

を直接検証し、finite-scale estimatorのbiasを調べる。

---

## B. Newton-degenerate capacity model

\[
F(x,y,z)
=
\frac{(x^2-1)^2}{4}
+
x^2(y^2+z^2)
+
y^6+y^2z^2+z^6.
\]

第一の直接予測は、

\[
\boxed{
\operatorname{cap}_T
\sim
C
T\log(1/T)e^{-1/(4T)}
}
\]

である。

exit-time predictionは、

\[
\boxed{
\mathbb E\tau
\sim
C'
\frac{T^{1/2}}{\log(1/T)}
e^{1/(4T)}.
}
\]

numerical testでは、

\[
\log\mathbb E\tau
-
\frac{1}{4T}
-
\frac12\log T
\]

を計算し、

\[
-\log\log\frac1T
+
c
\]

へ近づくかを見る。

---

## C. Branching and non-splittable models

第二段階以降では、

- monkey-saddle 型 branching gate、
- unstable mode と kernel singularity が分離できない analytic germ、

を構成する。

ここでは committor geometry 自体が新しい情報を持つかを検証する。

---

# 24. 経験的予言の優先順位

## Primary prediction

small-noise capacity prefactor。

第一モデルでは、

\[
\frac{
\operatorname{cap}_T
}{
T\log(1/T)e^{-1/(4T)}
}
\to
C.
\]

---

## Secondary prediction

exit-time negative logarithmic exponent。

\[
\log\tau
-
\frac{1}{4T}
-
\frac12\log T
\sim
-\log\log\frac1T
+
O(1).
\]

---

## Tertiary prediction

learning-model state transitionsのclock scaling。

LLC trajectoryとの時間的一致は、その後の観測仮説である。

本理論の勝負は、

> phase transitionが存在する

と言うことではない。

\[
\boxed{
\text{singular gate geometryからtransition prefactorを予測できるか}
}
\]

である。

---

# 25. 研究プログラム

## 第一撃

### A logarithmic Eyring–Kramers prefactor for a Newton-degenerate saddle

明示三次元potentialを扱う。

Newton-adapted transverse sectionを用いてsharp capacity lower boundを構成する。

目標：

\[
\operatorname{cap}_T
\sim
C
T\log(1/T)e^{-1/(4T)}.
\]

対応するexit-time lawを導く。

---

## 第二撃

### Split singular gates and transverse pole data

局所形

\[
F-H
=
-\frac{\mu s^2}{2}
+
K(v)
+
R(s,v)
\]

を考える。

\(K\) のresolution-of-singularities / Newton pole dataからcapacity asymptoticsを分類する条件を求める。

目標：

\[
\boxed{
\text{transverse pole data}
\longrightarrow
\text{capacity power/log data}
}
\]

を定理化する。

---

## 第三撃

### Branching and non-splittable singular capacities

fiberwise one-dimensional reductionが壊れるgateを扱う。

committor germ、harmonic splitting、gate topologyを含む新しいcapacity invariantを探索する。

これが完全なsingular capacity theoryの中心問題である。

---

## 第四撃

### Metastable Markov networks with singular weights and capacities

node weights

\[
Z_i(T)
\]

とedge capacities

\[
\operatorname{cap}_T(A_i,A_j)
\]

からreduced rate matrixを構成する。

Markov closure errorとtime-scale separationを評価する。

---

## 第五撃

### Learning-model bridge

reduced-rank regression、matrix factorization、deep linear networksなどで、algebraic strataとmetastable basinsの対応を調べる。

目標は、

\[
\boxed{
\text{learning stages}
\longleftrightarrow
\text{singular metastable transition network}
}
\]

を具体的モデルで検証することである。

---

## 第六撃

### Nonequilibrium singular geometry

state-dependent anisotropic diffusion

\[
dW_t
=
b(W_t)\,dt
+
\sqrt{\varepsilon}\sigma(W_t)\,dB_t
\]

に対し、sharp prefactor asymptoticsを研究する。

regular-variation exponentsが存在する条件を求める。

可逆極限でsingular mass/capacity calculusへ還元されるかを問う。

---

# 26. 本理論が既存理論と異なる点

## 特異学習理論との関係

SLTはsingular Laplace integralsをRLCT pole dataで分類する。

DSLTはそのmachineryを、まずwell massesの記述に用いる。

さらに本命として、

> capacity側に対応するsingular asymptotic classificationを構成できるか

を問う。

SLTを置換しない。

---

## Eyring–Kramers理論との関係

Morse、nonquadratic finite-type、Morse–Bott、convex-profileなどの退化Eyring–Kramers理論は既に広いclassを扱う。

したがって本理論は、

> 退化点で非Gaussian powerが出る

ことを新規性とは主張しない。

また、

\[
\lambda-\kappa
\]

という差構造そのものも新規性とは主張しない。

新規性候補は、

\[
\boxed{
\text{既存局所モデルの外にあるanalytic singular gatesについて、Newton/resolution dataをsharp capacity asymptoticsへ接続すること}
}
\]

である。

第一標的は、Newton-degenerate transverse geometryによるleading logarithmic capacity prefactorである。

---

## Developmental interpretabilityとの関係

LLC trajectoryはlocal singular complexityのprobeとなり得る。

DSLTはLLCをphase markerとして使うだけではない。

将来的には、

\[
\text{local geometry}
\to
\text{node mass / edge capacity data}
\to
\text{transition clock}
\]

という接続を目指す。

ただしLLCとtrue asymptotic pole dataの同一視はしない。

---

## Grokking理論との関係

grokkingをdelayed basin transitionとみなすこと自体は独自主張ではない。

DSLTの予言候補は、

\[
\frac{\text{barrier}}{\text{noise}}
+
\text{power exponent}\cdot\log(\text{noise})
+
\text{log multiplicity}\cdot\log\log(1/\text{noise})
\]

というclock scalingを、singular geometryから予測することである。

---

# 27. 正直に難しいところ

## 第一の障害

### Singular capacity

capacityはequilibrium potentialを含むPDE variational quantityである。

split gateではtransverse Laplace geometryへ還元できる可能性がある。

しかしbranching/non-splittable gateでは単純還元が壊れる。

これが最大の数学的障害である。

---

## 第二の障害

### Local versus global capacity data

capacity exponentがlocal gate germだけで決まる保証はない。

multiple gates、parallel channels、basin-pair geometryが同じexponential scaleで競合する可能性がある。

---

## 第三の障害

### Algebraic strata and metastable basins

algebraic stratificationがそのままdynamical state decompositionになるとは限らない。

両者の対応条件が必要である。

---

## 第四の障害

### Nonreversible SGD

RLCTがanisotropic SGD noiseの下で保存される保証はない。

stable power/log exponents自体が存在しない可能性もある。

---

## 第五の障害

### Finite-noise regime

理論は

\[
T\to0
\]

または

\[
\varepsilon\to0
\]

の漸近である。

deep learningのfinite-noise regimeではsubleading correctionが支配的な可能性がある。

---

## 第六の障害

### Local coefficient estimation

empirical LLC estimatorとtrue asymptotic pole dataの対応は自明ではない。

観測されたjumpを直ちにalgebraic transitionと呼んではならない。

---

# 28. 反証条件

## Falsification I

第一のNewton-degenerate gate modelで、

\[
T\log(1/T)e^{-1/(4T)}
\]

というpredicted capacity scalingが成立しない。

これは第一数学命題への直接反証である。

---

## Falsification II

split singular gate classにおいて、transverse singular Laplace pole dataとcapacity power/log dataの間に安定な対応が存在しない。

この場合、transverse pole classification programを棄却する。

---

## Falsification III

branching/non-splittable gatesでcapacity dataがlocal geometryから再現可能に分類できず、essentially global basin-pair informationへ依存する。

この場合、general local singular-capacity invariant構想を縮小する。

---

## Falsification IV

SGD-like diffusionでstable sub-Arrhenius power/log exponentsが存在しない。

この場合、nonequilibrium singular geometry構想を棄却する。

---

## Falsification V

learning-model stage transitionsでmetastable first-exit pictureとtime-scale separationが成立しない。

この場合、DSLTのlearning applicationは限定的small-noise regimeへ縮小される。

---

# 29. この理論が正しければ何が言えるか

第一に、過剰パラメータモデルの準安定状態占有を、

\[
L
\]

だけでも、

dimensionだけでも、

Hessian flatnessだけでもなく、

singular Gibbs mass dataとして記述できる。

第二に、transition clockを、

\[
\boxed{
\text{well mass}
\to
\text{capacity network}
\to
\text{sharp exit-time law}
}
\]

として定量化できる。

第三に、gate geometryがpower correctionだけでなく、

\[
\log(1/T)
\]

のmultiplicity correctionを持ち得るなら、学習時間に負のlogarithmic exponentのような新しいscaling signatureが現れる。

第四に、learning rate、batch size、noise injection、weight decayを、

\[
\boxed{
\text{quasipotential}
+
\text{state mass geometry}
+
\text{transition conductance geometry}
}
\]

を変形する操作として理解できる可能性がある。

第五に、developmental interpretabilityを、

> どの状態にいるか

の検出から、

> 次のtransition clockをどの幾何が支配しているか

の解析へ拡張できる可能性がある。

---

# 30. 最終命題

動的特異学習理論の最終目標は、

> なぜ学習はある表現に滞留し、なぜある時点で別の表現へ移り、なぜその順序と時間尺度が再現可能なのか

という問いに対し、

\[
\boxed{
\text{特異幾何}
+
\text{準安定性}
+
\text{確率的散逸}
}
\]

という一つの数学言語を与えることである。

そのため本理論では、

\[
\lambda_A-\kappa_{A,B}
\]

という metastable singularity gap を organizing notation として用いながら、数学的中心をより根本的な capacity classification problem に置く。

\[
\boxed{
\textbf{井戸側の特異幾何は質量を分類する。}
}
\]

\[
\boxed{
\textbf{ゲート側の特異幾何はcapacity、すなわち伝導度を分類する。}
}
\]

\[
\boxed{
\textbf{準安定時間は両者の比として現れる。}
}
\]

したがって、動的特異学習理論の最小の数学的問いは次である。

\[
\boxed{
\textbf{
SLT が singular Laplace integrals を RLCT pole data で分類するように、
singular transition capacities を分類する幾何学的不変量を構成できるか。
}
}
\]

第一標的は、

\[
\boxed{
\textbf{
Newton-degenerate transverse gate における leading logarithmic capacity prefactor
}
}
\]

である。

その先に、

\[
\boxed{
\textbf{
potential germ
+
committor geometry
+
gate topology
}
}
\]

を統合する singular capacity theory がある。
これが、動的特異学習理論の出発点である。