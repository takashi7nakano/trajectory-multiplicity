# TRAJECTORY_MULTIPLICITY_PROTOCOL v1.0 — FROZEN 2026-09-21

**Project:** Separate line from PHI_FLOW and from Paper I — *Does extensive-lesion robustness of a
feed-forward reactive process contain information not present in its one-time reactive statistics,
and if so, is the operative additional information route multiplicity (scalar) or route connectivity
(structural)?*
**Status:** v1.0, FROZEN 2026-09-21 after B's final audit of draft v0.3 (PASS FOR FREEZE; two
editorial corrections and one algorithm specification applied, no substantive change). Derived from
`9_Trajectory_Multiplicity/ASSESSMENT_2026-09-21_rev3.md` (conceptual review closed by B). Draft
history: v0.1 B MAJOR REVISION; v0.2 B MINOR REVISION BEFORE FREEZE; v0.3 B PASS FOR FREEZE. §12.
**Freeze rule.** No endpoint, threshold, family definition, predictor role, or kill criterion changes
after this freeze. Every value previously marked [draft] is now the registered value; the tag is kept
in the text as a record of which values were choices. Nothing has been computed at freeze.
Computation begins only on an explicit instruction from the PI, following §9. Format follows `8_PHI_FLOW_PROTOCOL_v1.0.md` and SHAPE_VARIABLES v2.0.1.
**Scope:** feed-forward (layered, no revisit) A→B Markov dynamics only. No claim about recurrent
dynamics, continuous-state systems, weak-noise limits, biological systems, or PHI_FLOW.
**Relation to PHI_FLOW:** none beyond the methodological form (registered hypotheses, kill criteria,
controls before data, B re-implementation). No quantity of PHI_FLOW appears here and no result is
transferred in either direction.
**Executes:** Pipeline A (all construction, predictors, lesion runs). Pipeline B: independent re-
implementation of the reactive kernels, all predictors, and the lesion response from the frozen
family specification, after A's results are frozen.

## 1. Registered question and hypotheses (from rev. 3 §1, unchanged)

Predictor classes by minimum reconstruction order of the unperturbed reactive process at resolution ℓ:
𝒬₁ = functionals of individual-time marginals only; 𝒬₂ = functionals requiring at least one joint
distribution across two times (itinerary quantities are 𝒬₂-reconstructible for a Markov reactive
process; noted, not tiered).

> **Q.** Across a family of feed-forward A→B chains with matched one-time reactive statistics, is the
> extensive-lesion collapse scale η_½(ℓ) predicted by the best 𝒬₁(ℓ) predictor? If not, which 𝒬₂(ℓ)
> predictor carries it — a structural connectivity statistic or a scalar (effective route number,
> coherence)?

  H₀     — η_½(ℓ) is 𝒬₁-predictable in the registered family (no variation beyond removal-set
           uncertainty resolved in the 𝒬₁-matched primary family, and no 𝒬₂ advantage in the
           secondary embedding). The line stops.
  H_red  — η_½(ℓ) is not 𝒬₁-predictable and is predicted by a cross-layer connectivity statistic of
           the reactive current at resolution ℓ; scalars are proxies. **Primary hypothesis.**
  H_mult — η_½(ℓ) is predicted by a scalar 𝒬₂ quantity better than by any connectivity statistic,
           over a family matched in 𝒬₁ and in cut statistics.

Written before computation (A and B, rev. 2–3): H_red expected. No novelty is claimed for "parallel
pathways produce robustness"; the only admissible claim is the outcome of the 𝒬₁/𝒬₂ separation.

## 2. Model class and construction rule (all fixed at freeze)

**Chain.** Layers t = 0..T, n positions per layer, x ∈ {0..n−1} (positions supply the coarse-graining
and the lesion footprint). Absorbing B after layer T; absorbing F reachable from every layer.
A = layer 0 with μ_A uniform over positions [draft]. Feed-forward: no revisits, so every single lesion
is exactly 𝒬₁-reducible (rev. 3 §4.3) and single lesions serve only as code checks.

**Construction rule (B11).** Specify the *reactive* transition kernels first, then embed:
1. Choose reactive marginals μ^R_t. Primary family: **uniform**, μ^R_t(x) = 1/n for all t [draft].
   (Any doubly stochastic kernel preserves it; committor-margin statistics are then identical across
   the family by construction, which is the intended matching, not a weakness — see §8.)
2. Choose reactive kernels P^R_t, doubly stochastic, from the menus of §3.
3. Embed by uniform leak: P_t(x,y) = (1−λ) P^R_t(x,y), P_t(x,F) = λ, for all t, x [draft]. Then the
   committor is q_t(x) = (1−λ)^{T−t} (constant within each layer), R₀ = (1−λ)^T, and the Doob-
   conditioned process is exactly P^R_t. R₀ is matched across the whole family identically.
   Secondary embedding (reported, not an endpoint) [draft]: position-dependent leak λ_t(x) with a
   registered profile, to verify that the endpoint ordering does not depend on the committor being
   flat within layers.

**Sizes** [draft]: n = 256, T = 64, λ such that R₀ = 0.5 (λ = 1 − 0.5^{1/64} ≈ 0.01077). Resolutions
ℓ ∈ {1, 2, 4, 8, 16} (blocks of ℓ consecutive positions; block index b = ⌊x/ℓ⌋). Nested dyadic
partition with origin at 0; partition-origin sensitivity check (C1) at shifts s ∈ {ℓ/4, ℓ/2, 3ℓ/4} for
ℓ ≥ 4.

## 3. The family (C2, C3, C5, C6) [all draft]

All kernels doubly stochastic on n = 256 positions with periodic boundary; a family member is one
kernel used at every layer (time-homogeneous) unless stated. Target family size N_fam ≥ 24 so that
rank correlations in §6 are meaningful; the menus below give ≥ 30.

  Loc(σ)      symmetric nearest-range walk: P^R(x,y) = 1/(2σ+1) for |x−y| ≤ σ;  σ ∈ {1, 2, 4, 8, 16, 32}
              (σ = 1: anchor System I, one broad local route; σ ≥ 16 approaches mean-field)
  Perm(π)     deterministic permutations: identity (anchor System II: n disconnected narrow routes),
              cyclic shift by s ∈ {1, 4, 16, 64}, and three fixed random permutations (seeds §7)
  Mix(c)      (1−c)·identity + c·uniform, c ∈ {0.02, 0.05, 0.1, 0.2, 0.5} (anchor System III at c = 0.1)
  Block(M,c)  M equal blocks of n/M positions; within-block uniform mixing with prob 1−c, uniform
              jump to another block with prob c; M ∈ {4, 16, 64}, c ∈ {0.02, 0.1}
  ShiftMix    (1−c)·shift(s) + c·uniform, s = 16, c ∈ {0.05, 0.2} (persistent transfer between blocks
              plus mixing)

**C2** is satisfied identically (uniform marginals, all members). **C3** (matched degree and flux
distributions): report the out-degree distribution and the edge-flux distribution of every member;
the analysis of §6 is repeated on the sub-family with matched ⟨k⟩ (Loc(σ) vs Block(M,c) at equal
support size) so that 1/⟨k⟩ alone cannot order the result. **C5** (matched cut structure, I vs III):
for the pair Loc(σ*) vs Mix(c*) choose (σ*, c*) [draft: σ* = 4, c* = 0.1] such that the cross-layer
min-cut at ℓ = 4 agrees within 5%; the residual Δη_½ is endpoint E4. **C6** (same marginals, same
scalar coherence, different wiring, cut width *not* matched): the pairs {identity, shift(16), random
permutation} (all C = T ln n), and {Mix(c), ShiftMix(c)} at equal c (equal C by construction since C
depends only on the row entropy of the kernel); per-block occupancy is uniform for all, so it is
matched. Endpoint E3.

**Controls before any family data (C4)** [draft]:
  K-a  Noise-inflated: Loc(1) with an extra independent coordinate z ∈ {0..15} resampled uniformly at
       every layer, lesions defined on x only. Every predictor computed on (x, z) must equal its value
       on x alone to 10⁻¹⁰ after itinerary classes are formed on x; raw path entropy differs by
       T ln 16 and must *not* be used.
  K-b  Timing control: Loc(1) with a gate at layer T/2 whose waiting time is bimodal (0 or 8 layers,
       equiprobable), realized by 8 extra waiting layers with identity kernel entered with prob 1/2.
       Itinerary N_eff at every ℓ must equal that of Loc(1) without the gate; fixed-time N_eff would
       double (reported to document the failure mode, not an endpoint).

## 4. Lesion family ℋ_ℓ (fixed)

Footprint ℓ: an element of ℋ_ℓ is a block index b ∈ {0..n/ℓ−1}. Extensive lesion at strength η: a
uniformly random subset of ⌈η n/ℓ⌉ block indices [draft].

**Time structure — primary: persistent** [draft]. A removed block is removed at *all* layers (a
spatial region lesioned for the whole duration). Secondary, reported: independent removal per
(layer, block) at the same η.

**Mode — primary: redirect** [draft]. For every layer t and every x, edges (x → y) with y in a removed
block are deleted and the row is renormalized proportionally over surviving successors (rev. 3 §4.3);
a row with no surviving successors becomes a dead end (→ F). Mass already in a removed block at layer
0 is lost. Secondary, reported: **dead-end** mode (a trajectory entering a removed block fails; no
renormalization). The difference R_redirect(η) − R_dead(η) is the **redirection gain** — the part of
robustness that the redirect intervention produces by moving mass onto surviving routes. It is
*intervention-dependent evidence for substitutability under proportional redirection*, not a direct
measure of substitutability (a different redirection rule gives a different gain); endpoint E5.

η grid: {0.05, 0.10, …, 0.95} [draft]; N_rem = 200 removal sets per (member, ℓ, η) [draft], each solved
exactly (backward committor recursion, O(T n²)). Collapse scale: η_½(ℓ) = the η at which
⟨R(η)⟩/R₀ crosses 1/2, linearly interpolated between grid points [draft]; also η_{1/e}. Error of η_½:
bootstrap over removal sets, 1000 resamples [draft].

## 5. Predictors (computed from the unperturbed reactive process, frozen before any lesion run)

All at resolution ℓ, on block indices, itinerary classes (sequence of blocks visited, consecutive
repeats removed) where itineraries are used.

  𝒬₁  P1  layerwise Hill numbers N^{(α)}_t of μ^R_t on blocks, α ∈ {0, 1, 2}: identical across the
          uniform family (= n/ℓ) — reported to show the matching, and non-trivial in the secondary
          embedding.
      P2  per-layer occupancy width W_t(ε): blocks with mass > ε/(n/ℓ), ε = 0.1 [draft]; min_t and mean.
      P3  committor-margin statistics: min over reactive paths of q (trivial in the primary
          embedding; non-trivial in the secondary).
  𝒬₂  P4  cross-layer reactive connectivity at resolution ℓ, **primary form = fractional max-flow**:
          on the layered block graph with edge capacity = reactive flux F^R_t(b, b′) and node capacity
          = block reactive mass μ^R_t(b), the A→B max-flow value normalized by total reactive mass
          (equivalently, by max-flow/min-cut, the smallest normalized cut of the reactive current
          separating A from B). Secondary form P4-int: the maximum number of node-disjoint A→B block-
          routes carrying flux > ε_f [draft ε_f = 10⁻⁴ of total]; reported, not used in E2.
      P5  adjacent-layer reactive edge network statistics: mean out-degree on blocks, and the
          persistence fraction (flux staying in the same block index across a layer).
      P6  scalar coherence C = Σ_t I(X_t; X_{t+1} | 𝒮) on blocks.
      P7  itinerary N_eff^{(1)}(ℓ) and N_eff^{(2)}(ℓ) (𝒬₂-reconstructible; computed by forward
          recursion on itinerary suffix classes, exact).
      P8  accessibility a(ℓ): for a single persistent block lesion in redirect mode, the fraction of
          the redirected mass that reaches B, averaged over blocks (a 𝒬₂ quantity by reconstruction
          order; it is the mean of the route-transition kernel's success column).
  Excluded as competitor: the perturbed committor; raw path entropy (see K-a).

Every predictor is checked for C1 invariance under block-label permutations (must be identical) and
under the partition-origin shifts of §2 (reported as a sensitivity; a predictor whose ranking of the
family changes by Spearman ρ < 0.9 across shifts is flagged and its endpoint contribution reported
with and without it) [draft].

## 6. Endpoints and decision rules [all draft; revised through v0.3]

**Pre-specified primary predictors (no best-of-class selection).** One predictor per role, fixed
before any lesion run:

  𝒬₁ primary:            P2  (per-layer occupancy width, min over layers)
  𝒬₂ structural primary: P4  (fractional max-flow of the reactive current)
  𝒬₂ scalar primary:     P7  (itinerary N_eff^{(1)}(ℓ) — the quantity the original hypothesis named)
  𝒬₂ scalar secondary:   P6  (scalar coherence C), reported alongside; see E2 rule.

All other predictors (P1, P3, P4-int, P5, P7^{(2)}, P8) are reported and enter only the exploratory
analysis E7. The v0.1 rule "largest |ρ| within class, evaluated on the same family" is withdrawn: with
unequal class sizes it gives 𝒬₂ a systematic selection advantage.

**Statistic.** Spearman ρ between η_½(ℓ) and the predictor across the family at each ℓ, with 90%
bootstrap CI over family members (resampling members) and over removal sets. Δ between two |ρ| values
has its CI from the same paired bootstrap.

**A structural fact that shapes E1.** In the primary family (C2 satisfied identically, uniform
marginals, uniform leak) *every* 𝒬₁ predictor is constant across members by construction — that is
what matching means. A correlation with a constant predictor is undefined, so E1 in the primary family
cannot be a ρ comparison; it is a variance test. The ρ comparison is meaningful only in the secondary
embedding, where the position-dependent leak makes μ^R_t, and hence P2, vary across the family. Both
are registered.

  E1a (H₀ in the primary family; F1′)  Does η_½(ℓ) vary across the family beyond removal-set
      uncertainty? Statistic [B, v0.2 review]: the excess variance

          V_excess(ℓ) = Var_members[η̂_½(ℓ)] − mean_members[s_m²(ℓ)],

      where η̂_½^{(m)}(ℓ) is member m's estimate over its N_rem removal sets and s_m²(ℓ) is that
      member's removal-set bootstrap variance. Uncertainty: hierarchical bootstrap — resample members
      with replacement, and within each resampled member resample its removal sets — 1000 outer
      resamples [draft]; 90% CI. Algorithm specification (B, final audit): when a member appears k
      times in an outer resample, each appearance receives its own independent inner removal-set
      resample; the k copies are treated as distinct members in computing Var_members and
      mean_members for that outer replicate. E1a *resolves* 𝒬₂-relevant variation at ℓ if the CI for V_excess(ℓ)
      lies above zero; otherwise no such variation is resolved at ℓ. F1′ fires at ℓ if no variation is
      resolved. H₀-consistent is declared if F1′ fires at every ℓ ∈ {2, 4, 8, 16}.
      Registered wording of that outcome: "no robustness variation beyond removal-set uncertainty was
      resolved in the registered 𝒬₁-matched family at the registered sizes" — not that one-time
      statistics fix extensive-lesion robustness in general. The line stops on that finding as a
      matter of allocation, not as a general claim.
  E1b (H₀ in the secondary embedding; F1′ ρ-form)  Δρ₁(ℓ) = |ρ(P4)| − |ρ(P2)|, both pre-specified.
      F1′ fires at ℓ if the 90% CI of Δρ₁(ℓ) includes 0 or Δρ₁ ≤ 0.10. Declared H₀-consistent if it
      fires at every ℓ. E1a and E1b are reported together; H₀-consistent is declared only if both
      hold, with the same restricted wording as E1a (a finding about the registered family, not a
      general statement); if they disagree, the disagreement is reported and no hypothesis is selected
      (K6).
  E2 (H_red vs H_mult; F1″)  Δρ₂(ℓ) = |ρ(P7)| − |ρ(P4)|, both pre-specified. F1″ fires (H_red) at ℓ if
      Δρ₂(ℓ) ≤ 0.10 or its 90% CI includes 0. H_red is declared if F1″ fires at every ℓ where E1 did
      not declare H₀. H_mult is declared only if Δρ₂ > 0.10 with a 90% CI excluding 0 at ≥ 2
      resolutions. P6 is not part of the decision: if P6 alone exceeds P4 by the same margin, this is
      reported as an exploratory finding for a future registered test, not as H_mult.
  E3 (C6 scalar-insufficiency; rewritten)  For each C6 pair (m, m′): D = η_½^{(m)}(ℓ) − η_½^{(m′)}(ℓ)
      with 90% CI from the removal-set bootstrap (paired over the same removal sets, which the shared
      block partition allows). A pair is "resolved" if the CI excludes 0. Scalar coherence is declared
      insufficient at ℓ if at least one pair is resolved with |D| > 3 σ_D; the number of resolved pairs
      and every D with its CI are reported. No family-range normalization and no half-pairs rule.
  E4 (C5 residual)  For the cut-matched pair, D with 90% CI at ℓ = 4. Expected: CI includes 0. A
      resolved |D| > 3σ_D is reported as "beyond cut statistics" and motivates a search for a further
      𝒬₂ or higher-order quantity; it does not by itself support H_mult.
  E5 (redirection gain; secondary)  G(η) = ⟨R_redirect − R_dead⟩/R₀ at η = η_½. Reported per member
      and correlated with P8. No decision rule; ordering by P8 is the expectation. Interpreted as
      evidence for substitutability under the registered redirection rule only.
  E6 (time-structure check; secondary)  η_½ under independent per-(layer, block) removal. Expected:
      Perm members become indistinguishable (survival (1−η)^T for all); reported to document why the
      persistent structure is primary.
  E7 (exploratory, out-of-sample; no decision weight)  The family is split into two halves by a
      registered seed (2026300; stratified by kernel menu). On half H₁ the predictor with largest |ρ|
      is selected within each class; its |ρ| is then evaluated on H₂, and the roles of H₁/H₂ are
      swapped. Reported as the only place where best-of-class appears; a predictor that wins out-of-
      sample but is not the pre-specified primary is a candidate for the next registered version, not
      a result.

Reporting rule: all ρ, all predictors, all members, all ℓ, both embeddings, both modes, both time
structures are reported whatever the endpoints show.

## 7. Code checks, seeds, frozen outputs [draft]

  T1  Reactive marginals: forward propagation of μ_A under P^R_t reproduces uniform to 10⁻¹².
  T2  Embedding: committor computed from the embedded chain equals (1−λ)^{T−t} to 10⁻¹²; Doob
      transform of the embedded chain reproduces P^R_t to 10⁻¹².
  T3  Entropy identity: H(Γ|𝒮) computed by the chain rule equals Σ_t H_t − C to 10⁻¹⁰ (rev. 3 §4.5).
  T4  Single-lesion exactness: ΔR from rev. 3 §4.3 formula equals the direct recomputation to 10⁻¹²
      for 100 random single-edge lesions per member (feed-forward exactness).
  T5  C1: every predictor identical under 10 random block-label permutations.
  T6  Controls K-a, K-b pass (§3).
Seeds: family random permutations 2026101–2026103; removal sets seed = 2026200000 + 10⁴·(member
index) + 10²·(ℓ index) + (η index), generator PCG64. Frozen outputs: `data/tm_family.json` (kernel
specifications), `data/tm_predictors_{ell}.csv` (frozen before lesion runs, SHA-256 in
`results/tm_manifest.md`), `data/tm_lesions_{member}_{ell}.npz`.

## 8. Kill criteria (any one → null report, no reinterpretation)

  K1  Any of T1–T6 fails → pipeline invalid; fix, bump protocol version, rerun checks before any
      family data are examined.
  K2  C2 matching violated for any member (marginal deviation > 10⁻¹⁰) → member excluded.
  K3  F1′ fires at every ℓ (E1a and E1b) → H₀-consistent in the registered family; the line stops;
      report written with the restricted wording of E1a.
  K4  F1″ fires → H_mult dead; H_red reported (this is the expected outcome and is not a failure of
      the study).
  K5  A predictor fails C1 (T5) → that predictor is excluded; the study continues.
  K6  E1a and E1b disagree (one declares H₀-consistent, the other not) at any ℓ, or the primary and
      secondary embeddings give opposite signs of Δρ₂ at any ℓ → the result is declared embedding-
      dependent and no hypothesis is selected.

## 9. Execution order

Freeze code and seeds → T1–T5 → controls K-a, K-b (T6) → construct family, compute and freeze all
predictors at all ℓ (manifest) → lesion runs (primary time structure and mode first) → E1 → E2 → E3,
E4 → E5, E6 → secondary embedding, secondary time structure and mode → Pipeline A results frozen →
Pipeline B re-implementation from `data/tm_family.json`; PASS iff every predictor agrees to 10⁻⁸ and
every η_½ agrees within bootstrap error.

## 10. What is not claimed

Anything about recurrent dynamics (a cyclic toy is required first, rev. 3 §1 Scope); anything about
continuous-state, chaotic, or weak-noise systems; anything about biological systems, brains, identity,
or consciousness; that "parallel pathways produce robustness" (established, §2.13 of rev. 3); that
route multiplicity is a new concept — the only claim available is which 𝒬 class, and within 𝒬₂
which kind of quantity, predicts extensive-lesion robustness after 𝒬₁ is matched; anything about
PHI_FLOW or Paper I.

## 11. Draft choices for B (collected)

  D1  Uniform reactive marginals as the primary family (committor flat within layers by construction;
      secondary embedding as the check). Alternative: a non-uniform marginal family from the start.
  D2  Persistent (spatial) block lesions as primary; independent per-(layer, block) as secondary.
  D3  [accepted by B with change, applied] Redirect primary, dead-end secondary; "redirection gain"
      E5 as intervention-dependent evidence for substitutability.
  D4  Collapse scale η_½ by interpolation (no sharp threshold at finite T); η_{1/e} alongside.
  D5  [rewritten after B] Pre-specified primary predictors P2 / P4 / P7 (P6 secondary scalar, no
      decision weight); E1 split into E1a (variance test in the C2-matched primary family, where 𝒬₁
      is constant by construction) and E1b (ρ-form in the secondary embedding); best-of-class only in
      the out-of-sample exploratory E7. [Settled by B on v0.2: P7 remains primary scalar, P6
      secondary; E1a accepted conceptually with the excess-variance statistic and hierarchical
      bootstrap of v0.3.]
  D6  [accepted by B with change, applied] Fractional max-flow primary; integer route count
      secondary (P4-int), not used in E2.
  D7  Sizes n = 256, T = 64, R₀ = 0.5, ℓ up to 16, N_rem = 200.
  D8  [accepted] K6 as a registered kill; restated in v0.2 to cover E1a/E1b disagreement.
  D9  [accepted by B] E3 decision: "at least one C6 pair resolved with |D| > 3σ_D", all D and CIs
      reported.
  D10 [new, v0.3] Hierarchical bootstrap for V_excess: 1000 outer (member) resamples with inner
      removal-set resampling; 90% CI. Only the resample count is a draft value.

## 12. Change log

- v0.1 DRAFT (2026-09-21): first draft from assessment rev. 3 after B's closure of the conceptual
  review. Not registered. Nothing computed.
- v0.2 DRAFT (2026-09-21): after B's review of v0.1. (1) "substitution gain" → "redirection gain",
  interpreted as intervention-dependent evidence (D3). (2) P4 primary form = fractional max-flow;
  integer route count secondary (D6). (3) Best-of-class selection withdrawn from all decision rules;
  primary predictors pre-specified as P2 / P4 / P7, P6 secondary scalar without decision weight;
  best-of-class confined to the out-of-sample exploratory E7 (D5, blocking issue). (4) E1 split into
  E1a (variance test) and E1b (ρ-form), because in the C2-matched primary family every 𝒬₁ predictor is
  constant by construction and a ρ comparison is undefined — a latent tautology in v0.1 that A caught
  while implementing B's change. (5) E3 threshold replaced by within-pair difference with bootstrap CI
  (D9). (6) K6 restated. Not registered. Nothing computed.
- v0.3 DRAFT (2026-09-21): after B's review of v0.2. (1) E1a statistic replaced: variance ratio V vs 1
  → excess variance V_excess = Var_members(η̂_½) − mean_members(s_m²), hierarchical bootstrap over
  members and removal sets, resolution iff the 90% CI lies above zero (D10). (2) Registered wording of
  a null E1a/E1b outcome restricted to "no robustness variation beyond removal-set uncertainty was
  resolved in the registered 𝒬₁-matched family", propagated to H₀ and K3. (3) §11 items settled by B
  marked as such. No other change. Not registered. Nothing computed.
- v1.0 FROZEN (2026-09-21): B final audit of v0.3 PASS FOR FREEZE. Editorial only: K1 version-
  specific wording removed; §6 heading updated; hierarchical-bootstrap rule for repeated members
  specified (algorithm specification, not an endpoint change). Frozen. Nothing computed.
