# Lean Formalization Design for Dynamical Singular Learning Theory

## 1. Purpose

This document specifies the formalization strategy for the Lean component of the Dynamical Singular Learning Theory (DSLT) research program.

The purpose of the Lean project is not to formalize the entirety of singular learning theory, metastability theory, stochastic differential equations, or Eyring–Kramers theory from first principles.

The initial purpose is narrower:

> **to identify, isolate, and mechanically verify the logical and asymptotic core of DSLT, while making every unformalized mathematical dependency explicit.**

The central candidate principle of DSLT is that the sub-Arrhenius power appearing in a metastable transition time is governed by a difference between a well-side singular Laplace exponent and a transition-capacity exponent.

Symbolically,

[
\delta_{A\to B}
:=
\lambda_A-\kappa_{A,B}.
]

The proposed transition-time asymptotic is

[
\mathbb E\tau_{A\to B}
\sim
C
T^{\lambda_A-\kappa_{A,B}}
\left(\log\frac1T\right)^{m_A-r_{A,B}}
\exp\left(
\frac{H_{A,B}-L_A}{T}
\right).
]

The Lean formalization must distinguish two fundamentally different tasks:

1. proving that the above asymptotic follows logically from suitable well-mass, capacity, and mass–capacity hypotheses;
2. proving that those hypotheses actually hold for a concrete metastable diffusion.

The first task belongs to the initial formalization program.

The second task is a substantially deeper mathematical problem and is initially treated as an explicit collection of proof obligations.

---

## 2. Research Philosophy

The formalization follows the principle

[
\text{claim}
\longrightarrow
\text{dependency decomposition}
\longrightarrow
\text{formal proof obligations}
\longrightarrow
\text{kernel-checked consequences}.
]

Lean is not used as a substitute for mathematical assumptions.

In particular, the formalization must never hide an unresolved mathematical problem behind a custom axiom.

The role of Lean is to answer questions such as:

* What exactly is assumed?
* What exactly follows?
* Which positivity conditions are required?
* Which logarithms are well-defined?
* Where is an asymptotic equivalence used?
* Is a sign or exponent incorrect?
* Is an allegedly new theorem only an algebraic consequence of known hypotheses?
* Which mathematical statements remain genuinely unproved?

The project therefore treats formalization as a tool for **proof-obligation discovery**, not merely proof certification.

---

## 3. Central Formal Structure

The initial DSLT bridge has three mathematical inputs.

### 3.1 Well-mass asymptotic

For a positive temperature parameter (T\to0^+), suppose a well mass (M(T)) satisfies

[
M(T)
\sim
C_M
e^{-L/T}
T^\lambda
\left(\log\frac1T\right)^{m-1}.
]

Here:

* (C_M>0) is a leading constant;
* (L) is the well energy level;
* (\lambda) is the well singular Laplace exponent;
* (m) is the logarithmic multiplicity parameter.

In a concrete DSLT application, (M(T)) is expected to arise from a local Gibbs mass,

[
M(T)
====

\int_A e^{-L(w)/T},dw.
]

The general derivation of this asymptotic from RLCT or resolution-of-singularities theory is initially outside the formalization scope.

---

### 3.2 Capacity asymptotic

Suppose a transition capacity (K(T)) satisfies

[
K(T)
\sim
C_K
e^{-H/T}
T^\kappa
\left(\log\frac1T\right)^{r-1}.
]

Here:

* (C_K>0) is a leading capacity constant;
* (H) is a communication height;
* (\kappa) is a capacity asymptotic exponent;
* (r) is a logarithmic multiplicity parameter.

In a concrete reversible Langevin model, (K(T)) is expected to be a potential-theoretic capacity.

The existence, locality, and geometric classification of (\kappa) are not assumed to be solved by DSLT.

The general theory of singular capacity asymptotics is a separate research target.

---

### 3.3 Mass–capacity exit relation

Suppose an exit-time quantity (E(T)) satisfies

[
E(T)
\sim
\frac{M(T)}{K(T)}.
]

In metastability theory, the precise form of this statement depends on the process, the metastable sets, and the choice of initial law.

Therefore the Lean project must not silently identify an arbitrary expected hitting time with a mass–capacity quotient.

The mass–capacity relation is initially represented as an explicit hypothesis.

Its concrete justification must later be linked to a precise metastability theorem with clearly stated assumptions.

---

### 3.4 Abstract DSLT bridge

From the preceding hypotheses, DSLT predicts

[
E(T)
\sim
\frac{C_M}{C_K}
e^{(H-L)/T}
T^{\lambda-\kappa}
\left(\log\frac1T\right)^{m-r}.
]

The metastable singularity gap is defined by

[
\delta:=\lambda-\kappa.
]

The first formal objective is to prove the abstract implication

[
\boxed{
\text{well asymptotic}
+
\text{capacity asymptotic}
+
\text{mass–capacity relation}
\Longrightarrow
\text{singular exit-time asymptotic}.
}
]

This result is referred to in the Lean project as the **Abstract Bridge Theorem**.

The Abstract Bridge Theorem is not itself a proof of the full DSLT metastability conjecture.

It is a formally checked reduction of that conjecture to explicit analytic and probabilistic proof obligations.

---

## 4. Formalization Boundaries

### 4.1 Initially in scope

The first formalization stage includes:

1. elementary asymptotic equivalence lemmas;
2. quotient laws for asymptotic equivalents;
3. exponential-factor cancellation;
4. power-exponent subtraction;
5. logarithmic-multiplicity subtraction;
6. the logarithmic clock law;
7. the local free-energy expansion;
8. the coarse-grained free-energy identity;
9. the definition of the metastable singularity gap;
10. the Morse cancellation check;
11. exact quartic-well scaling;
12. exact separable monomial-well scaling;
13. the Abstract Bridge Theorem.

---

### 4.2 Initially out of scope

The following topics are deliberately excluded from the first formalization stage:

* general resolution of singularities;
* Hironaka-style desingularization;
* general RLCT theory;
* meromorphic continuation of local zeta functions;
* general singular Laplace asymptotics;
* construction of overdamped Langevin diffusions;
* general stochastic differential equation theory;
* quasi-stationary distributions;
* general potential-theoretic metastability;
* Dirichlet-form capacity theory;
* sharp Eyring–Kramers theorems;
* Witten Laplacians;
* singular transition gates;
* singular capacity localization;
* Freidlin–Wentzell theory;
* nonreversible diffusion;
* SGD diffusion approximation;
* grokking models.

These exclusions are methodological rather than conceptual.

An out-of-scope topic may later become a formalization phase if it becomes necessary for a DSLT theorem.

---

## 5. Core Definitions

### 5.1 Temperature domain

All small-temperature asymptotics are taken as

[
T\to0^+.
]

The Lean representation should use a filter describing convergence to zero from the positive side.

No theorem involving

[
\log T,
\qquad
T^\alpha,
\qquad
\log\log\frac1T
]

may rely on an implicit assumption that (T>0).

Positivity and eventual domain conditions must be made explicit.

---

### 5.2 Asymptotic equivalence

The preferred mathematical notation is

[
f(T)\sim g(T)
\quad
(T\to0^+).
]

The Lean implementation should use Mathlib's asymptotic-equivalence infrastructure whenever practical.

A project-local notation or wrapper may be introduced only if it improves readability without changing the underlying semantics.

The project should avoid defining an independent, incompatible notion of asymptotic equivalence.

---

### 5.3 Singular asymptotic profile

A reusable abstraction should represent functions of the form

[
C
e^{-E/T}
T^\alpha
\left(\log\frac1T\right)^q.
]

A possible conceptual structure is

```lean
structure SingularProfile where
  constant : ℝ
  energy : ℝ
  power : ℝ
  logPower : ℝ
```

This exact Lean representation is not mandatory.

Before introducing such a structure, the implementation should determine whether direct functions and lemmas produce simpler proofs.

The formalization should prefer mathematical transparency over premature abstraction.

---

### 5.4 Metastable singularity gap

The fundamental DSLT definition is

[
\delta
======

\lambda-\kappa.
]

Conceptually:

* (\lambda) is a well-side exponent;
* (\kappa) is a transition-capacity exponent;
* (\delta) is associated with a directed transition.

The formalization should preserve the distinction between state quantities and transition quantities.

A future graph-level formalization may encode

[
\delta_{i\to j}
===============

\lambda_i-\kappa_{ij}.
]

The initial implementation may use real numbers without introducing basin or graph types.

---

## 6. Initial Theorem Inventory

### DSLT-T001: Asymptotic Quotient Law

Given

[
f\sim f_0,
\qquad
g\sim g_0,
]

under suitable eventual nonvanishing assumptions,

[
\frac{f}{g}
\sim
\frac{f_0}{g_0}.
]

This theorem should preferably reuse existing Mathlib asymptotic lemmas rather than reprove general facts unnecessarily.

The DSLT-specific theorem should then specialize this law to singular asymptotic profiles.

---

### DSLT-T002: Singular Profile Quotient

Assume

[
M(T)
\sim
C_Me^{-L/T}T^\lambda
\left(\log\frac1T\right)^{m-1},
]

and

[
K(T)
\sim
C_Ke^{-H/T}T^\kappa
\left(\log\frac1T\right)^{r-1}.
]

Then

[
\frac{M(T)}{K(T)}
\sim
\frac{C_M}{C_K}
e^{(H-L)/T}
T^{\lambda-\kappa}
\left(\log\frac1T\right)^{m-r}.
]

The formal proof must explicitly verify the exponent identities

[
\lambda-\kappa
]

and

[
(m-1)-(r-1)=m-r.
]

This theorem is the algebraic asymptotic core of the metastable singularity-gap law.

---

### DSLT-T003: Logarithmic Clock Law

Suppose

[
\tau(T)
\sim
C
T^\delta
\left(\log\frac1T\right)^q
e^{\Delta/T},
]

with suitable positivity assumptions.

Then

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
o(1).
]

The formalization must carefully distinguish between:

* exact logarithmic identities for the asymptotic reference profile;
* transfer of asymptotic equivalence through the logarithm;
* the proof that the ratio converges to (1);
* eventual positivity required for taking logarithms.

The theorem must not treat

[
f\sim g
]

as if it directly implied

[
\log f\sim\log g.
]

The desired result concerns an additive (o(1)) difference, not multiplicative asymptotic equivalence of logarithms.

---

### DSLT-T004: Local Free-Energy Expansion

Let

[
Z(T)
\sim
C
e^{-L/T}
T^\lambda
\left(\log\frac1T\right)^{m-1}.
]

Define

[
F(T):=-T\log Z(T).
]

Then

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
o(T).
]

This theorem is an important regression test because the signs of the logarithmic-multiplicity and constant terms are easy to state incorrectly.

The expected signs are

[
-(m-1)T\log\log\frac1T
]

and

[
-T\log C.
]

---

### DSLT-T005: Coarse-Grained Free-Energy Identity

Let (I) be a finite state space.

For positive weights (Z_i), define

[
F_i=-T\log Z_i,
]

[
Z_{\mathrm{cg}}=\sum_iZ_i,
]

and

[
\pi_i=\frac{Z_i}{Z_{\mathrm{cg}}}.
]

For a probability vector (p), define

[
G(p,T)
======

\sum_i p_iF_i
+
T\sum_i p_i\log p_i.
]

Then

[
\boxed{
G(p,T)
======

## T D_{\mathrm{KL}}(p\Vert\pi)

T\log Z_{\mathrm{cg}}.
}
]

The finite-state theorem should explicitly handle conventions involving (0\log0).

The initial implementation may impose strict positivity on (p_i) if this substantially simplifies Phase 1.

A later theorem may remove that restriction.

---

### DSLT-T006: Morse Cancellation

Define

[
\delta=\lambda-\kappa.
]

If

[
\lambda=\frac d2
]

and

[
\kappa=\frac d2,
]

then

[
\delta=0.
]

Consequently,

[
T^\delta=1
]

for (T>0).

This theorem is mathematically elementary but serves as a formal regression test for the claim that the DSLT power correction disappears in the regular Morse limit.

---

### DSLT-T007: Quartic-Well Scaling

For (T>0), define

[
Z_4(T)
======

\int_{\mathbb R}
e^{-x^4/T},dx.
]

Prove the exact identity

[
\boxed{
Z_4(T)
======

T^{1/4}Z_4(1).
}
]

This theorem should be proved by a rigorous change of variables.

The result is exact and should not be weakened to an asymptotic statement.

The formalization must also establish the required integrability.

---

### DSLT-T008: Separable Monomial-Well Scaling

Let

[
L(x)
====

\sum_{i=1}^d
c_ix_i^{2k_i},
]

where

[
c_i>0
]

and

[
k_i\in\mathbb N_{>0}.
]

Define

[
Z(T)
====

\int_{\mathbb R^d}
\exp\left(
-\frac{L(x)}{T}
\right),dx.
]

Prove

[
Z(T)
====

C
T^{\sum_i1/(2k_i)}
]

for a positive constant (C) independent of (T).

The theorem should preferably derive the result from coordinate rescaling and product integration.

The associated well exponent is

[
\lambda
=======

\sum_i\frac1{2k_i}.
]

This theorem is the first mechanically verified example in which a non-Gaussian local geometry produces a non-Morse temperature exponent.

---

### DSLT-T009: Abstract Bridge Theorem

Assume:

1. `WellMassAsymptotic`;
2. `CapacityAsymptotic`;
3. `MassCapacityExitRelation`.

Then prove:

[
E(T)
\sim
\frac{C_M}{C_K}
e^{(H-L)/T}
T^{\lambda-\kappa}
\left(\log\frac1T\right)^{m-r}.
]

Equivalently, with

[
\delta=\lambda-\kappa,
]

prove

[
E(T)
\sim
C
T^\delta
\left(\log\frac1T\right)^{m-r}
e^{(H-L)/T}.
]

The theorem must make all three mathematical inputs visible in its theorem signature or in named hypotheses.

The project must not encode any of these inputs as a global axiom.

---

## 7. Proof-Dependency Graph

The initial dependency structure is:

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

Independent formal branches are:

```text
Laplace Scaling
    ├── Quartic Well
    └── Separable Monomial Well
```

and

```text
Free Energy
    ├── Local Free-Energy Expansion
    └── Coarse-Grained Free-Energy Identity
```

The intended long-term dependency graph is:

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

Only the central implication is initially formalized.

The three incoming mathematical theories remain independent proof obligations.

---

## 8. Claim Status Classification

Every mathematically substantive DSLT claim should be assigned one of the following statuses.

### `FORMALIZED`

A Lean theorem exists.

The theorem contains no `sorryAx`.

No custom mathematical axioms are required.

---

### `EXTERNAL`

The claim is taken from existing mathematical literature.

The repository must record:

* an exact bibliographic reference;
* the theorem or proposition number where possible;
* the hypotheses of the cited result;
* the mapping between the source notation and DSLT notation.

An `EXTERNAL` claim is not Lean-verified merely because it is cited.

---

### `DERIVED`

The claim is a mathematical consequence of `FORMALIZED` and/or explicitly stated `EXTERNAL` claims.

If the derivation is central to DSLT, it should eventually become `FORMALIZED`.

---

### `CONJECTURE`

The claim is mathematically unproved within the project and has not been identified as a known theorem.

Conjectures must not be represented as Lean axioms.

---

### `EMPIRICAL`

The claim concerns numerical experiments, machine-learning observations, or parameter estimation.

Empirical claims are outside Lean proof status.

---

### `OPEN`

The correct mathematical statement is not yet known.

This status is preferable to prematurely writing a conjecture when the expected theorem form itself is uncertain.

---

## 9. Axiom and `sorry` Policy

The formalization adopts a strict axiom policy.

### 9.1 No custom mathematical axioms

The project must not introduce declarations such as

```lean
axiom singularCapacityAsymptotic : ...
```

for unresolved DSLT claims.

Doing so would allow downstream proofs to be kernel-checked while hiding the central mathematical problem.

---

### 9.2 `sorry` is allowed only during active development

Temporary use of `sorry` is permitted in local work-in-progress branches.

A theorem containing or depending on `sorryAx` must not be:

* listed as formally proved;
* cited in a paper as Lean-verified;
* included in a formal-verification result table as complete.

---

### 9.3 Public theorem audit

Public DSLT theorems should be audited using

```lean
#print axioms DSLT.theoremName
```

where practical.

The output should be reviewed before a theorem is marked `FORMALIZED`.

---

### 9.4 Trusted foundations

Dependencies on Lean's standard foundational axioms and Mathlib's accepted classical infrastructure must be documented honestly.

The project does not claim constructive mathematics unless a theorem is specifically proved constructively.

The objective is to prevent DSLT-specific assumptions from being disguised as foundational dependencies.

---

## 10. File Organization

The initial Lean source tree is:

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

The directory names represent mathematical responsibilities.

A file should not be placed in `Metastability/` merely because DSLT is motivated by metastability.

For example, the quotient of two asymptotic profiles belongs in `Asymptotics/`.

---

## 11. Naming Conventions

Lean declaration names should describe mathematical content rather than the history of the project.

Preferred:

```lean
asymptotic_singularProfile_div
localFreeEnergy_expansion
coarseGrainedFreeEnergy_eq_kl
singularityGap_eq_zero_of_morse
quarticPartition_scaling
monomialPartition_scaling
exit_asymptotic_of_mass_capacity
```

Avoid names such as:

```lean
mainTheorem
importantLemma
newDSLTResult
proof1
claimA
```

The project namespace should be

```lean
namespace DSLT
```

with suitable nested namespaces when helpful.

---

## 12. Treatment of External Theorems

A central methodological problem is connecting informal literature results to Lean statements.

The project must not write

> by Eyring–Kramers

inside the conceptual argument without recording the actual hypotheses.

For every external theorem used in DSLT, create an entry in

```text
literature/claims-map.md
```

containing:

1. DSLT claim identifier;
2. mathematical statement;
3. status;
4. source;
5. source theorem number;
6. source assumptions;
7. DSLT assumptions;
8. assumption-mapping notes;
9. known gaps.

Example:

```markdown
## DSLT-E001: Morse-gate capacity asymptotic

Status: EXTERNAL / not yet formally imported

DSLT form:

cap_T(A, B) ~ C exp(-H/T) T^(d/2)

Required mapping questions:

- Is the capacity normalized with or without the Gibbs partition function?
- Does the source use epsilon or T?
- Is the Dirichlet form multiplied by T?
- Is the saddle unique?
- Are competing gates excluded?
- What regularity is required of L?
- Is the result local or global?
```

The assumption mapping is part of the research result.

---

## 13. Numerical Experiments and Lean

Lean does not validate asymptotic claims by numerical fitting.

Numerical experiments serve a different purpose.

For example, a quartic-well experiment may numerically verify that

[
\frac{Z_4(T)}{T^{1/4}}
]

is approximately constant.

Lean should prove the exact scaling law.

The numerical experiment tests implementation and intuition.

The formal theorem tests deduction.

The two must not be conflated.

---

## 14. AI-Assisted Formalization Policy

The project explicitly permits AI systems to assist with:

* Lean code generation;
* theorem decomposition;
* Mathlib search;
* error interpretation;
* proof repair;
* informal-to-formal translation;
* literature-query generation;
* mathematical criticism.

However, AI-generated Lean code has no privileged epistemic status.

A proof is considered formally verified only if it is accepted by Lean under the project's axiom policy.

AI-generated mathematical explanations may contain errors even when associated Lean code compiles.

Conversely, compiling code may prove a theorem whose assumptions are too strong or mathematically irrelevant.

Therefore every theorem requires two reviews:

1. **formal review:** does Lean check the statement and proof?
2. **semantic review:** is the statement actually the intended DSLT claim?

This distinction is mandatory.

---

## 15. Expert-Feedback Loop

The project is designed around the iterative cycle

[
\text{AI hypothesis generation}
\to
\text{formal decomposition}
\to
\text{Lean verification}
\to
\text{expert criticism}
\to
\text{claim revision}
\to
\text{reformalization}.
]

Expert feedback should be recorded as research input rather than silently incorporated.

When a substantive criticism changes:

* a definition;
* a theorem hypothesis;
* a claimed novelty;
* a proof dependency;
* the interpretation of an exponent,

the change should be documented in the research log.

Where appropriate, a design decision record should explain the revision.

The objective is to preserve the development history of the theory.

---

## 16. Initial Milestones

### Milestone 1: Formal asymptotic core

Complete without `sorry`:

* DSLT-T001;
* DSLT-T002;
* DSLT-T003;
* DSLT-T004;
* DSLT-T005;
* DSLT-T006.

Success means the algebraic and logarithmic structure of the DSLT formulas has been mechanically checked.

---

### Milestone 2: Exact singular-well toy models

Complete:

* DSLT-T007;
* DSLT-T008.

Success means at least one non-Morse well exponent has been derived from an exact, kernel-checked integral scaling argument.

---

### Milestone 3: Abstract Bridge Theorem

Complete DSLT-T009.

Success means Lean verifies that the metastable singularity-gap asymptotic follows from three explicitly named mathematical inputs.

At this stage, the repository must clearly state:

> The full DSLT metastability theorem has not yet been proved.

The remaining proof obligations must be visible.

---

### Milestone 4: External theorem audit

Identify precise literature results corresponding to:

* singular well-mass asymptotics;
* Morse-gate capacity asymptotics;
* mass–capacity exit-time relations.

For each result, perform a hypothesis and normalization audit.

The main research question at this stage is:

> Do existing theorems already imply the first DSLT theorem, or is an additional bridge/localization argument required?

---

### Milestone 5: First concrete metastable theorem

State one concrete theorem for a restricted model class, preferably:

[
\text{degenerate analytic well}
+
\text{unique Morse index-one gate}.
]

The theorem should use an explicitly defined process, basin, initial law, capacity convention, and asymptotic regime.

Only after this theorem is proved mathematically should the project claim a theorem-level result in dynamical singular learning theory.

---

## 17. Definition of Phase-1 Success

Phase 1 succeeds if the repository establishes the following without custom DSLT axioms or `sorry` dependencies:

1. the DSLT singular-profile quotient formula is correct;
2. the power exponent is exactly (\lambda-\kappa);
3. the logarithmic exponent is exactly (m-r);
4. the logarithmic clock law is correctly derived;
5. the local free-energy expansion has the correct signs;
6. the coarse-grained free-energy identity has the correct normalization sign;
7. the Morse case gives zero singularity gap;
8. quartic and separable monomial wells exhibit formally verified non-Gaussian scaling;
9. the Abstract Bridge Theorem is formally proved from explicit hypotheses.

Phase 1 does **not** establish that the hypotheses of the Abstract Bridge Theorem hold for general Langevin metastability.

Its purpose is to transform the DSLT proposal from an informal narrative into a precise dependency graph of mathematical claims.

---

## 18. Definition of Theoretical Success

The first theorem-level success of DSLT requires more than the Abstract Bridge Theorem.

A concrete theorem must establish, for a mathematically specified stochastic system, that:

[
\mathbb E\tau_{A\to B}
\sim
C
T^{\lambda_A-d/2}
\left(\log\frac1T\right)^{m_A-1}
e^{(L(\sigma)-L_A)/T}
]

under explicit assumptions involving a singular well and a Morse transition gate.

The project must determine whether this statement is:

1. already an immediate consequence of existing results;
2. a new bridge theorem obtained by combining existing theories;
3. false without additional assumptions;
4. a genuinely new metastability theorem.

The result of this classification is itself a major objective of the project.

---

## 19. Long-Term Formalization Questions

The following questions are intentionally left open.

### Q1. Can RLCT data be represented naturally in Lean?

Can the relevant pole data be formalized without first formalizing a large part of resolution-of-singularities theory?

---

### Q2. What is the correct formal object for a capacity exponent?

Is

[
\kappa_{A,B}
]

a property of:

* a local potential germ;
* a potential and committor pair;
* a Dirichlet-form germ;
* a basin pair;
* a gate network?

The formal type should not prejudge this research question.

---

### Q3. Can a singular capacity profile be defined by regular variation?

Rather than assuming a complete power-log expansion, should the first invariant be defined by the existence of

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

\kappa?
]

This may provide a weaker and more stable formal target.

---

### Q4. Are singularity gaps constrained by reversibility?

For a reversible metastable Markov reduction, detailed balance may impose consistency relations among:

[
\lambda_i,
\qquad
\kappa_{ij},
\qquad
\delta_{i\to j}.
]

Possible cycle or edge-reversal identities should be investigated.

---

### Q5. What survives in nonreversible dynamics?

For small-noise nonreversible diffusions, does a stable exponent

[
\gamma_{A\to B}
]

exist?

Can it be decomposed as

[
\gamma_{A\to B}
===============

\chi_A-\chi^{\mathrm{gate}}_{A,B}?
]

No such decomposition is assumed in the initial formalization.

---

## 20. Final Principle

The Lean project should make it impossible to confuse the following statements:

> The DSLT formula is algebraically consistent.

> The DSLT formula follows from explicit asymptotic hypotheses.

> Existing mathematical literature proves those hypotheses.

> The first DSLT theorem has been proved.

> Singular capacity theory has been constructed.

> DSLT applies to SGD.

These are different levels of mathematical achievement.

The formalization is successful when those levels are represented by different declarations, different claim statuses, and different proof dependencies.

The guiding principle of the project is:

[
\boxed{
\text{Do not formalize confidence. Formalize assumptions and consequences.}
}
]

The role of Lean in DSLT is therefore not to certify a grand theory prematurely.

Its role is to expose exactly where the grand theory begins.