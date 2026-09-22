# TRAJECTORY_MULTIPLICITY_PROTOCOL v1.1 — FREEZE RECORD

**凍結日:** 2026-09-21
**対象:** `PROTOCOL_v1.1.md`
**SHA-256:** `4a34a4e1dea12adbb94bb5c067e1691a34d4301632f4600974ab890d2660a9dd`
**状態:** **FROZEN. 凍結時点で outcome は何も計算していない。**

## レビュー履歴（すべて lesion outcome 前）

| 段階 | B の裁定 |
|---|---|
| v1.0 | FROZEN → P4 退化発見 → **SUPERSEDED BEFORE LESION OUTCOME COMPUTATION** |
| 解析評価 | **ANALYTIC ASSESSMENT PASS** |
| REPAIR SPEC v1.1 | **MINOR REVISION**（closure 対称化・C3′・C6/E3 を control 留置） |
| v1.1 DRAFT | **MINOR REVISION BEFORE FREEZE**（E2 三帰結・K4・E1 降格・T8 注記） |
| v1.1 RC rev1 | 上記 4 点を反映 |
| v1.1 RC rev2 | **ONE BLOCKING CLARIFICATION** に対応：member → family 集約規則（existential rule）を登録、K3/K4/K6 更新 |
| v1.1 | **PASS FOR FREEZE**（2026-09-21） |

B PASS FOR FREEZE 全文（2026-09-21）:
「The existential member-to-family rule is logically aligned with the registered universal
sufficiency question. The disclosure in §6.2.6 is appropriate prospective transparency, not
outcome-based rule selection, because it uses only already-computed moments and analytic
properties of the registered closure while D_{k,j}(ℓ) remains uncomputed. No further substantive
protocol changes are requested.」

## 凍結された設計の要点

- **primary intervention:** dead-end（persistent block lesion）。redirect は secondary。
- **registered question:** dead-end robustness curve は E[V(ℓ)] だけで再現されるか、
  P(V(ℓ)) の高次情報を要するか。情報階層 E[V] ⊂ {E[V], Var(V)} ⊂ P(V)。
- **primary estimator:** 除去集合を解析的に積分した軌道推定量（removal-set noise ゼロ）、M = 10⁶。
- **oracle:** 厳密 P(V) が得られる 13 member（置換 8 直接列挙、Mix 5 占有数閉形式）。implementation validation。
- **validation:** removal-set simulation（v1.0 seed 体系、N_rem = 200）。
- **primary endpoint E2:** member ごとの sup ノルム曲線距離 D_{k,j}(ℓ)、τ_R = 10⁻²。
  member → family は existential rule（family H₀ ⟺ 全 m_j=0；H_dist ⟺ ∃j m_j≥2；さもなくば inconclusive）。
- **撤回:** P4（族内定数）、P7（全族厳密計算不能）、E5（gain 解析的既知）、C5。
- **降格:** E1（記述的診断、decision weight なし）。C6/E3（strong control）。
- **新しい structural graph statistic は探さない**（endpoint が P(V) のみで決まることが証明済み）。

## 凍結宣言

registered removal seed は未消費。`data/tm_lesions_*.npz`・`data/tm_curves_*.npz`・
`data/tm_PV_*.npz` はいずれも不存在。族の η_½・robustness curve・距離 D_{k,j} は一つも計算していない。

**計算は §9 の実行順序に従い、PI（中野）の明示的な指示があってはじめて開始する。**
Pipeline A 結果凍結後に Pipeline B が `data/tm_family.json` から独立再実装する。

## 検証済みの事項（凍結時点）

- T1–T5: 27/27 member 通過（`tm_checks_v11.json`、最大残差 2.1e-12）。
- T7: 解析恒等式 4 本を族外トイ系で総当たり検証（最大 1.2e-14）。
- T9: closure-1/closure-2 が全 27 member × 5 解像度で well-posed（失敗 0）。
- feasibility: 厳密 P(V) の 3 経路（列挙・占有数・標本）が互いに、かつ厳密 E[V] と一致。

## 未実施（凍結後、PI 指示で行う）

- T6（対照 K-a, K-b の実 checkでの検証）
- T8 の curve 照合（oracle 厳密 P(V) vs 推定量、グリッド点ごと 4 SE）
- E[V]・Var(V) の全 member 凍結 manifest
- 軌道標本 M = 10⁶ 抽出、P̂(V) 凍結
- E2 → E3 → E6 → E7 → E1、secondary、removal-set validation
