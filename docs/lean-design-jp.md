# 動的特異学習理論の Lean 形式化設計

## 1. 目的

本文書は、動的特異学習理論（Dynamical Singular Learning Theory; DSLT）研究プログラムにおける Lean 形式化の方針を定めるものである。

本 Lean プロジェクトの目的は、特異学習理論、準安定性理論、確率微分方程式、Eyring–Kramers 理論の全体を第一原理から形式化することではない。

初期段階の目的は、より限定的である。

> **DSLT の論理的・漸近的な数学的核を分離し、機械的に検証すると同時に、未形式化・未証明の数学的依存関係をすべて明示する。**

DSLT の中心的な原理候補は、準安定遷移時間に現れる sub-Arrhenius な温度 power が、井戸側の特異 Laplace 指数と transition capacity 指数の差によって支配される、というものである。

象徴的には、

[
\delta_{A\to B}
:=
\lambda_A-\kappa_{A,B}.
]

提案されている遷移時間漸近は、

[
\mathbb E\tau_{A\to B}
\sim
C
T^{\lambda_A-\kappa_{A,B}}
\left(\log\frac1T\right)^{m_A-r_{A,B}}
\exp\left(
\frac{H_{A,B}-L_A}{T}
\right)
]

である。

Lean 形式化では、次の二つの課題を厳密に区別しなければならない。

1. 適切な well mass、capacity、mass–capacity relation の仮定から、上記漸近式が論理的に導かれることを証明する。
2. 具体的な準安定拡散過程に対して、それらの仮定が実際に成立することを証明する。

第一の課題は、初期形式化プログラムの対象とする。

第二の課題は、はるかに深い数学的問題であり、初期段階では明示的な proof obligation、すなわち「今後証明すべき命題群」として扱う。

---

## 2. 研究哲学

本形式化では、

[
\text{主張}
\longrightarrow
\text{依存関係の分解}
\longrightarrow
\text{形式的 proof obligation}
\longrightarrow
\text{kernel により検査された帰結}
]

という原則を採用する。

Lean を数学的仮定の代用品として使用してはならない。

特に、未解決の数学問題を独自 axiom の背後へ隠してはならない。

Lean の役割は、例えば次の問いに答えることである。

* 正確には何を仮定しているのか。
* その仮定から正確には何が導かれるのか。
* どの positivity 条件が必要なのか。
* どの logarithm が well-defined である必要があるのか。
* asymptotic equivalence はどこで使用されているのか。
* 符号や指数に誤りはないか。
* 新定理と思われているものが、実際には既知の仮定からの代数的帰結にすぎないのではないか。
* 真に未証明の数学的主張はどこに残っているのか。

したがって本プロジェクトでは、形式化を単なる proof certification の手段としてではなく、**proof obligation を発見するための道具**として扱う。

---

## 3. 中心的な形式構造

初期 DSLT bridge は、三つの数学的入力を持つ。

### 3.1 Well mass の漸近

正の温度パラメータ

[
T\to0^+
]

を考える。

well mass (M(T)) が、

[
M(T)
\sim
C_M
e^{-L/T}
T^\lambda
\left(\log\frac1T\right)^{m-1}
]

を満たすと仮定する。

ここで、

* (C_M>0) は leading constant、
* (L) は well の energy level、
* (\lambda) は well singular Laplace exponent、
* (m) は logarithmic multiplicity parameter

である。

具体的な DSLT の応用では、

[
M(T)
====

\int_A e^{-L(w)/T},dw
]

という局所 Gibbs mass から (M(T)) が生じることを想定する。

RLCT や resolution of singularities から、この漸近を一般的に導出する理論は、初期形式化の対象外とする。

---

### 3.2 Capacity の漸近

transition capacity (K(T)) が、

[
K(T)
\sim
C_K
e^{-H/T}
T^\kappa
\left(\log\frac1T\right)^{r-1}
]

を満たすと仮定する。

ここで、

* (C_K>0) は leading capacity constant、
* (H) は communication height、
* (\kappa) は capacity asymptotic exponent、
* (r) は logarithmic multiplicity parameter

である。

具体的な可逆 Langevin 系では、(K(T)) は potential-theoretic capacity に対応することを想定する。

しかし DSLT は、

* (\kappa) が存在すること、
* (\kappa) が局所的な量であること、
* (\kappa) が特異幾何学的に分類可能であること

を、既に解決済みとは仮定しない。

一般的な singular capacity asymptotics の構成は、独立した研究課題である。

---

### 3.3 Mass–capacity exit relation

exit-time quantity (E(T)) が、

[
E(T)
\sim
\frac{M(T)}{K(T)}
]

を満たすと仮定する。

準安定性理論において、この命題の正確な形は、

* 確率過程、
* metastable set の定義、
* 初期分布

に依存する。

したがって Lean プロジェクトは、任意の expected hitting time を暗黙に mass–capacity quotient と同一視してはならない。

mass–capacity relation は、初期段階では明示的な仮定として表現する。

具体的な正当化は、後に明確な仮定を持つ既存の準安定性定理、あるいは新しい定理へ接続しなければならない。

---

### 3.4 抽象 DSLT bridge

以上の仮定から DSLT は、

[
E(T)
\sim
\frac{C_M}{C_K}
e^{(H-L)/T}
T^{\lambda-\kappa}
\left(\log\frac1T\right)^{m-r}
]

を予測する。

metastable singularity gap を、

[
\delta:=\lambda-\kappa
]

と定義する。

最初の形式化目標は、

[
\boxed{
\text{well asymptotic}
+
\text{capacity asymptotic}
+
\text{mass–capacity relation}
\Longrightarrow
\text{singular exit-time asymptotic}
}
]

という抽象的含意を証明することである。

Lean プロジェクトでは、この結果を **Abstract Bridge Theorem** と呼ぶ。

Abstract Bridge Theorem 自体は、完全な DSLT metastability conjecture の証明ではない。

これは DSLT の中心予想を、明示された解析的・確率論的 proof obligation へ還元する、形式的に検査された reduction theorem である。

---

## 4. 形式化の境界

### 4.1 初期段階で対象とするもの

第一形式化段階では、以下を対象とする。

1. 基本的な asymptotic equivalence lemma
2. asymptotic equivalent の quotient law
3. exponential factor の cancellation
4. power exponent の subtraction
5. logarithmic multiplicity の subtraction
6. logarithmic clock law
7. local free-energy expansion
8. coarse-grained free-energy identity
9. metastable singularity gap の定義
10. Morse cancellation の確認
11. quartic well の厳密な scaling
12. separable monomial well の厳密な scaling
13. Abstract Bridge Theorem

---

### 4.2 初期段階では対象外とするもの

以下は意図的に第一形式化段階から除外する。

* 一般的な resolution of singularities
* Hironaka 型 desingularization
* 一般 RLCT 理論
* local zeta function の meromorphic continuation
* 一般 singular Laplace asymptotics
* overdamped Langevin diffusion の構成
* 一般確率微分方程式理論
* quasi-stationary distribution
* 一般 potential-theoretic metastability
* Dirichlet-form capacity theory
* sharp Eyring–Kramers theorem
* Witten Laplacian
* singular transition gate
* singular capacity localization
* Freidlin–Wentzell theory
* nonreversible diffusion
* SGD diffusion approximation
* grokking model

これらを除外する理由は概念的なものではなく、方法論的なものである。

対象外の問題であっても、DSLT の定理に必要となれば、将来の形式化 phase として追加する。

---

## 5. 中心定義

### 5.1 温度領域

すべての small-temperature asymptotics は、

[
T\to0^+
]

として扱う。

Lean では、正側から零へ収束する filter を用いるべきである。

次の量を含む定理では、

[
\log T,
\qquad
T^\alpha,
\qquad
\log\log\frac1T
]

暗黙に (T>0) を仮定してはならない。

positivity 条件、および eventual domain condition は明示する。

---

### 5.2 漸近同値

数学的な基本記法として、

[
f(T)\sim g(T)
\qquad
(T\to0^+)
]

を用いる。

Lean 実装では、可能な限り Mathlib の asymptotic equivalence infrastructure を使用する。

project-local な notation や wrapper は、underlying semantics を変更せず、可読性を実質的に改善する場合にのみ導入する。

独自の非互換な asymptotic equivalence 定義を作ることは避ける。

---

### 5.3 Singular asymptotic profile

次の形の関数を再利用可能な抽象として表現することを検討する。

[
C
e^{-E/T}
T^\alpha
\left(\log\frac1T\right)^q.
]

概念的には、例えば次の structure が考えられる。

```lean
structure SingularProfile where
  constant : ℝ
  energy : ℝ
  power : ℝ
  logPower : ℝ
```

ただし、この Lean 表現を必須とはしない。

このような structure を導入する前に、直接的な関数定義と lemma の方が proof を単純にできないかを確認する。

本形式化では、premature abstraction より数学的透明性を優先する。

---

### 5.4 Metastable singularity gap

DSLT の基本定義は、

[
\delta
======

\lambda-\kappa
]

である。

概念的には、

* (\lambda) は well-side exponent、
* (\kappa) は transition-capacity exponent、
* (\delta) は directed transition に付随する量

である。

形式化では、state quantity と transition quantity の区別を保存しなければならない。

将来 graph-level formalization を行う場合、

[
\delta_{i\to j}
===============

\lambda_i-\kappa_{ij}
]

と表現できる。

初期実装では basin type や graph type を導入せず、実数上の定義としてよい。

---

## 6. 初期定理一覧

### DSLT-T001: Asymptotic Quotient Law

[
f\sim f_0,
\qquad
g\sim g_0
]

を仮定する。

適切な eventual nonvanishing 条件の下で、

[
\frac{f}{g}
\sim
\frac{f_0}{g_0}
]

を証明する。

可能な限り、一般的事実を独自に再証明せず、Mathlib の既存 asymptotic lemma を使用する。

その後 DSLT 固有の theorem として、singular asymptotic profile に特殊化する。

---

### DSLT-T002: Singular Profile Quotient

次を仮定する。

[
M(T)
\sim
C_Me^{-L/T}T^\lambda
\left(\log\frac1T\right)^{m-1},
]

[
K(T)
\sim
C_Ke^{-H/T}T^\kappa
\left(\log\frac1T\right)^{r-1}.
]

このとき、

[
\frac{M(T)}{K(T)}
\sim
\frac{C_M}{C_K}
e^{(H-L)/T}
T^{\lambda-\kappa}
\left(\log\frac1T\right)^{m-r}
]

を証明する。

formal proof では、

[
\lambda-\kappa
]

および、

[
(m-1)-(r-1)=m-r
]

という exponent identity を明示的に確認する。

この theorem は metastable singularity-gap law の代数的・漸近的核である。

---

### DSLT-T003: Logarithmic Clock Law

次を仮定する。

[
\tau(T)
\sim
C
T^\delta
\left(\log\frac1T\right)^q
e^{\Delta/T}.
]

適切な positivity 条件の下で、

[
\log\tau(T)
===========

\frac{\Delta}{T}
+
\delta\log T
+
q\log\log\frac1T
+
\log C
+
o(1)
]

を証明する。

形式化では、次を慎重に区別する。

* asymptotic reference profile に対する exact logarithmic identity
* logarithm を通した asymptotic information の transfer
* ratio が (1) へ収束すること
* logarithm を取るための eventual positivity

特に、

[
f\sim g
]

から直接、

[
\log f\sim\log g
]

を導いてはならない。

目的の結論は logarithm の multiplicative asymptotic equivalence ではなく、additive (o(1)) difference である。

---

### DSLT-T004: Local Free-Energy Expansion

次を仮定する。

[
Z(T)
\sim
C
e^{-L/T}
T^\lambda
\left(\log\frac1T\right)^{m-1}.
]

[
F(T):=-T\log Z(T)
]

と定義する。

このとき、

[
F(T)
====

L
+
\lambda T\log\frac1T
--------------------

## (m-1)T\log\log\frac1T

T\log C
+
o(T)
]

を証明する。

この theorem は重要な regression test である。

logarithmic multiplicity 項、および constant 項の符号は誤記しやすい。

正しい符号は、

[
-(m-1)T\log\log\frac1T
]

および、

[
-T\log C
]

である。

---

### DSLT-T005: Coarse-Grained Free-Energy Identity

有限状態空間 (I) を考える。

正の重み (Z_i) に対し、

[
F_i=-T\log Z_i,
]

[
Z_{\mathrm{cg}}=\sum_iZ_i,
]

[
\pi_i=\frac{Z_i}{Z_{\mathrm{cg}}}
]

と定義する。

probability vector (p) に対して、

[
G(p,T)
======

\sum_i p_iF_i
+
T\sum_i p_i\log p_i
]

と定義する。

このとき、

[
\boxed{
G(p,T)
======

## T D_{\mathrm{KL}}(p\Vert\pi)

T\log Z_{\mathrm{cg}}
}
]

を証明する。

有限状態 theorem では、

[
0\log0
]

の convention を明示的に処理する必要がある。

Phase 1 の形式化を大幅に単純化できる場合、初期実装では

[
p_i>0
]

を仮定してもよい。

その制約は後の theorem で除去する。

---

### DSLT-T006: Morse Cancellation

[
\delta=\lambda-\kappa
]

と定義する。

もし、

[
\lambda=\frac d2
]

かつ、

[
\kappa=\frac d2
]

なら、

[
\delta=0
]

を証明する。

従って (T>0) に対し、

[
T^\delta=1.
]

数学的には初等的な theorem である。

しかし DSLT の power correction が regular Morse limit で消失するという主張に対する formal regression test として使用する。

---

### DSLT-T007: Quartic-Well Scaling

(T>0) に対して、

[
Z_4(T)
======

\int_{\mathbb R}
e^{-x^4/T},dx
]

と定義する。

厳密な等式、

[
\boxed{
Z_4(T)
======

T^{1/4}Z_4(1)
}
]

を証明する。

proof は rigorous change of variables による。

これは exact result であり、asymptotic statement へ弱めてはならない。

必要な integrability も形式化する。

---

### DSLT-T008: Separable Monomial-Well Scaling

[
L(x)
====

\sum_{i=1}^d
c_ix_i^{2k_i}
]

とする。

ここで、

[
c_i>0
]

かつ、

[
k_i\in\mathbb N_{>0}.
]

[
Z(T)
====

\int_{\mathbb R^d}
\exp\left(
-\frac{L(x)}{T}
\right),dx
]

と定義する。

正の定数 (C) が存在して、

[
Z(T)
====

C
T^{\sum_i1/(2k_i)}
]

を満たすことを証明する。

(C) は (T) に依存しない。

可能であれば、coordinate rescaling と product integration により証明する。

対応する well exponent は、

[
\lambda
=======

\sum_i\frac1{2k_i}.
]

これは non-Gaussian local geometry が non-Morse temperature exponent を生成することを機械的に検証する最初の例である。

---

### DSLT-T009: Abstract Bridge Theorem

以下を仮定する。

1. `WellMassAsymptotic`
2. `CapacityAsymptotic`
3. `MassCapacityExitRelation`

このとき、

[
E(T)
\sim
\frac{C_M}{C_K}
e^{(H-L)/T}
T^{\lambda-\kappa}
\left(\log\frac1T\right)^{m-r}
]

を証明する。

また、

[
\delta=\lambda-\kappa
]

と定義すれば、

[
E(T)
\sim
C
T^\delta
\left(\log\frac1T\right)^{m-r}
e^{(H-L)/T}
]

と書ける。

この theorem の signature、または named hypothesis の中に、三つの数学的入力をすべて可視化しなければならない。

いずれの入力も global axiom として encode してはならない。

---

## 7. Proof Dependency Graph

初期 dependency structure は次の通りである。

```text
Asymptotic Quotient Law
          │
          ▼
Singular Profile Quotient
          │
          ├───────────────┐
          ▼               ▼
Abstract Bridge      Log Clock Law
          │
          ▼
Metastable Singularity-Gap Consequence
```

独立した形式化 branch として、

```text
Laplace Scaling
    ├── Quartic Well
    └── Separable Monomial Well
```

および、

```text
Free Energy
    ├── Local Free-Energy Expansion
    └── Coarse-Grained Free-Energy Identity
```

を持つ。

長期的には、

```text
Singular Laplace Theory
          │
          ▼
Well Mass Asymptotic ─────────┐
                              │
Capacity Theory               │
          │                   │
          ▼                   ▼
Capacity Asymptotic ──► Abstract Bridge
                              ▲
                              │
Potential-Theoretic           │
Metastability ────────────────┘
```

という dependency graph を想定する。

初期段階では中央の implication のみを形式化する。

三つの入力数学理論は、独立した proof obligation として残す。

---

## 8. Claim Status の分類

数学的に重要な DSLT claim には、以下の status のいずれかを付与する。

### `FORMALIZED`

Lean theorem が存在する。

theorem は `sorryAx` を含まない。

独自の数学的 axiom を必要としない。

---

### `EXTERNAL`

既存数学文献から採用した claim。

repository には以下を記録する。

* 正確な bibliographic reference
* 可能であれば theorem または proposition number
* cited result の hypotheses
* source notation と DSLT notation の対応

単に citation が存在するだけでは、`EXTERNAL` claim は Lean verified とはみなさない。

---

### `DERIVED`

`FORMALIZED` claim、および明示された `EXTERNAL` claim から数学的に導かれる claim。

DSLT の中心的 derivation であれば、最終的には `FORMALIZED` へ移す。

---

### `CONJECTURE`

project 内では証明されておらず、既知定理としても確認されていない claim。

conjecture を Lean axiom として表現してはならない。

---

### `EMPIRICAL`

数値実験、機械学習上の観測、parameter estimation に関する claim。

Lean proof status の対象外である。

---

### `OPEN`

正しい数学的 statement 自体がまだ分かっていない問題。

期待される theorem form が不明確な段階で、早すぎる conjecture を作るより `OPEN` を使用する。

---

## 9. Axiom および `sorry` policy

本形式化では厳格な axiom policy を採用する。

### 9.1 独自の数学的 axiom を禁止する

次のような declaration を、未解決 DSLT claim のために導入してはならない。

```lean
axiom singularCapacityAsymptotic : ...
```

これを行うと、下流 theorem は kernel-checked されても、中心数学問題が隠蔽される。

---

### 9.2 `sorry` は active development 中のみ許可する

local work-in-progress branch では、一時的な `sorry` の使用を許可する。

ただし `sorryAx` を含む、または間接的に依存する theorem を、

* formal proof 完了として一覧へ載せる
* 論文で Lean-verified と記述する
* formal-verification result table で complete とする

ことは禁止する。

---

### 9.3 Public theorem audit

公開する DSLT theorem は、可能な限り、

```lean
#print axioms DSLT.theoremName
```

により audit する。

theorem を `FORMALIZED` とする前に出力を確認する。

---

### 9.4 Trusted foundations

Lean の標準的 foundational axiom、および Mathlib が採用している classical infrastructure への依存は正直に記録する。

特定 theorem を constructive に証明していない限り、本 project は constructive mathematics を主張しない。

目的は DSLT 固有の assumption が foundational dependency を装って隠れることを防ぐことである。

---

## 10. ファイル構成

初期 Lean source tree は次の通りとする。

```text
lean/
├── lakefile.toml
├── lean-toolchain
└── DSLT/
    ├── Asymptotics/
    │   ├── Basic.lean
    │   ├── Quotient.lean
    │   └── LogClockLaw.lean
    │
    ├── Definitions/
    │   └── SingularityGap.lean
    │
    ├── FreeEnergy/
    │   ├── LocalExpansion.lean
    │   └── CoarseGrained.lean
    │
    ├── Laplace/
    │   ├── QuarticWell.lean
    │   └── MonomialWell.lean
    │
    ├── Metastability/
    │   └── AbstractBridge.lean
    │
    └── DSLT.lean
```

directory name は数学的責任範囲を表す。

DSLT が metastability を動機としているという理由だけで、任意の file を `Metastability/` へ置いてはならない。

例えば、二つの asymptotic profile の quotient に関する theorem は `Asymptotics/` に属する。

---

## 11. Naming Convention

Lean declaration name は、project の歴史ではなく数学的内容を記述する。

推奨例：

```lean
asymptotic_singularProfile_div
localFreeEnergy_expansion
coarseGrainedFreeEnergy_eq_kl
singularityGap_eq_zero_of_morse
quarticPartition_scaling
monomialPartition_scaling
exit_asymptotic_of_mass_capacity
```

避けるべき例：

```lean
mainTheorem
importantLemma
newDSLTResult
proof1
claimA
```

project namespace は、

```lean
namespace DSLT
```

とする。

必要に応じて nested namespace を使用する。

---

## 12. External Theorem の扱い

DSLT における重要な方法論的問題の一つは、informal literature result と Lean statement を接続することである。

本 project は、

> by Eyring–Kramers

のような記述だけで conceptual argument を進めてはならない。

使用する external theorem ごとに、

```text
literature/claims-map.md
```

へ以下を記録する。

1. DSLT claim identifier
2. 数学的 statement
3. status
4. source
5. source theorem number
6. source assumptions
7. DSLT assumptions
8. assumption-mapping note
9. known gap

例：

```markdown
## DSLT-E001: Morse-gate capacity asymptotic

Status: EXTERNAL / not yet formally imported

DSLT form:

cap_T(A, B) ~ C exp(-H/T) T^(d/2)

Required mapping questions:

- capacity は Gibbs partition function で正規化されているか。
- source は epsilon と T のどちらを使用しているか。
- Dirichlet form に T が掛かっているか。
- saddle は unique か。
- competing gate は除外されているか。
- L にはどの regularity が要求されるか。
- result は local theorem か global theorem か。
```

assumption mapping 自体を研究成果の一部として扱う。

---

## 13. 数値実験と Lean

Lean は numerical fitting によって asymptotic claim を検証するものではない。

数値実験は異なる役割を持つ。

例えば quartic well experiment では、

[
\frac{Z_4(T)}{T^{1/4}}
]

が数値的にほぼ一定になることを確認できる。

Lean は exact scaling law を証明する。

数値実験は implementation と intuition を検査する。

formal theorem は deduction を検査する。

両者を混同してはならない。

---

## 14. AI-assisted formalization policy

本 project では AI system を以下の用途に使用することを明示的に許可する。

* Lean code generation
* theorem decomposition
* Mathlib search
* error interpretation
* proof repair
* informal-to-formal translation
* literature query generation
* mathematical criticism

ただし、AI が生成した Lean code に特権的な epistemic status はない。

proof が formally verified とされるのは、project の axiom policy の下で Lean が受理した場合のみである。

AI-generated mathematical explanation は、対応する Lean code が compile していても誤りを含み得る。

逆に compile する code が、過度に強い assumptions や数学的に無意味な statement を証明している可能性もある。

したがって、各 theorem には二種類の review が必要である。

1. **formal review:** Lean が statement と proof を検査するか。
2. **semantic review:** statement は本当に意図された DSLT claim か。

この区別は必須である。

---

## 15. Expert-feedback loop

本 project は、次の iterative cycle を前提として設計する。

[
\text{AIによる仮説生成}
\to
\text{形式的分解}
\to
\text{Lean検証}
\to
\text{専門家による批判}
\to
\text{claim修正}
\to
\text{再形式化}
]

専門家からの feedback は、黙って project へ吸収するのではなく、research input として記録する。

重要な criticism により、

* definition
* theorem hypothesis
* novelty claim
* proof dependency
* exponent interpretation

のいずれかが変更された場合、その変更を research log に記録する。

必要に応じて design decision record を作成し、revision の理由を説明する。

目的は、理論の development history を保存することである。

---

## 16. 初期 milestone

### Milestone 1: Formal asymptotic core

以下を `sorry` なしで完成させる。

* DSLT-T001
* DSLT-T002
* DSLT-T003
* DSLT-T004
* DSLT-T005
* DSLT-T006

成功条件は、DSLT formula の代数的・logarithmic structure が機械的に検査されたことである。

---

### Milestone 2: Exact singular-well toy models

以下を完成させる。

* DSLT-T007
* DSLT-T008

成功条件は、少なくとも一つの non-Morse well exponent が、exact kernel-checked integral scaling argument から導かれることである。

---

### Milestone 3: Abstract Bridge Theorem

DSLT-T009 を完成させる。

成功条件は、metastable singularity-gap asymptotic が、明示的に名前を付けられた三つの数学的入力から従うことを Lean が検証することである。

この段階では repository に明確に次を記述する。

> DSLT の完全な metastability theorem はまだ証明されていない。

残る proof obligation を可視化する。

---

### Milestone 4: External theorem audit

以下に対応する正確な既存文献上の theorem を特定する。

* singular well-mass asymptotics
* Morse-gate capacity asymptotics
* mass–capacity exit-time relation

各 result について hypothesis と normalization を audit する。

この段階の中心研究問題は、

> 既存定理だけで第一 DSLT theorem は既に導かれるのか。それとも追加の bridge theorem または localization argument が必要なのか。

である。

---

### Milestone 5: 最初の具体的 metastable theorem

限定された model class に対して、具体的 theorem を一つ定式化する。

第一候補は、

[
\text{degenerate analytic well}
+
\text{unique Morse index-one gate}
]

である。

theorem では、以下を明示的に定義する。

* stochastic process
* basin
* initial law
* capacity convention
* asymptotic regime

この theorem が数学的に証明された後にのみ、project は dynamical singular learning theory における theorem-level result を主張する。

---

## 17. Phase 1 の成功条件

以下を、独自 DSLT axiom および `sorry` dependency なしで確立した場合、Phase 1 は成功とする。

1. DSLT singular-profile quotient formula が正しい。
2. power exponent が正確に (\lambda-\kappa) である。
3. logarithmic exponent が正確に (m-r) である。
4. logarithmic clock law が正しく導出される。
5. local free-energy expansion の符号が正しい。
6. coarse-grained free-energy identity の normalization sign が正しい。
7. Morse case で singularity gap が zero になる。
8. quartic well と separable monomial well が kernel-checked non-Gaussian scaling を示す。
9. Abstract Bridge Theorem が explicit hypotheses から形式的に証明される。

Phase 1 は、Abstract Bridge Theorem の assumptions が一般 Langevin metastability に対して成立することを示すものではない。

その目的は DSLT proposal を informal narrative から、明確な dependency graph を持つ数学的 claim 群へ変換することである。

---

## 18. 理論的成功の定義

DSLT の最初の theorem-level success には、Abstract Bridge Theorem 以上のものが必要である。

数学的に明示された stochastic system に対して、

[
\mathbb E\tau_{A\to B}
\sim
C
T^{\lambda_A-d/2}
\left(\log\frac1T\right)^{m_A-1}
e^{(L(\sigma)-L_A)/T}
]

を、singular well と Morse transition gate に関する explicit assumptions の下で確立しなければならない。

project はこの statement が、

1. 既存結果から直ちに従うのか。
2. 既存理論を接続する新しい bridge theorem なのか。
3. 追加 assumptions なしでは偽なのか。
4. genuinely new metastability theorem なのか。

を分類する。

この分類結果自体を project の主要目的の一つとする。

---

## 19. 長期的な形式化上の問い

以下の問題は意図的に open とする。

### Q1. RLCT data を Lean で自然に表現できるか

resolution-of-singularities theory の大部分を先に形式化せず、必要な pole data を形式化できるか。

---

### Q2. Capacity exponent の正しい formal object は何か

[
\kappa_{A,B}
]

は何の性質なのか。

* local potential germ
* potential と committor の pair
* Dirichlet-form germ
* basin pair
* gate network

のいずれなのか。

formal type を設計する段階で、この研究問題の答えを先取りしてはならない。

---

### Q3. Singular capacity profile を regular variation によって定義できるか

最初から完全な power-log expansion を仮定する代わりに、

[
\lim_{T\to0^+}
\frac{
\log\left(
e^{H/T}\operatorname{cap}_T(A,B)
\right)
}{
\log T
}
=

\kappa
]

の存在によって最初の invariant を定義するべきか。

これは、より弱く安定した formal target になる可能性がある。

---

### Q4. Reversibility は singularity gap に制約を与えるか

reversible metastable Markov reduction では detailed balance により、

[
\lambda_i,
\qquad
\kappa_{ij},
\qquad
\delta_{i\to j}
]

の間に consistency relation が生じる可能性がある。

cycle identity や edge-reversal identity の存在を調査する。

---

### Q5. Nonreversible dynamics では何が残るか

small-noise nonreversible diffusion において、安定した exponent

[
\gamma_{A\to B}
]

は存在するか。

また、

[
\gamma_{A\to B}
===============

\chi_A-\chi^{\mathrm{gate}}_{A,B}
]

という分解は可能か。

初期形式化では、このような decomposition を一切仮定しない。

---

## 20. 最終原則

Lean project は、次の statement が混同されることを防がなければならない。

> DSLT formula は代数的に整合している。

> DSLT formula は明示された漸近仮定から従う。

> 既存数学文献がその仮定を証明している。

> 第一 DSLT theorem が証明された。

> singular capacity theory が構成された。

> DSLT が SGD に適用できる。

これらは、数学的達成度の異なる段階である。

本形式化が成功したと言えるのは、それぞれが、

* 異なる declaration
* 異なる claim status
* 異なる proof dependency

として表現されている場合である。

本 project の guiding principle は次である。

[
\boxed{
\textbf{確信を形式化するな。仮定と帰結を形式化せよ。}
}
]

Lean の DSLT における役割は、壮大な理論を早すぎる段階で認証することではない。

その役割は、

> **壮大な理論が、数学として本当に始まる地点を露出すること**

である。
