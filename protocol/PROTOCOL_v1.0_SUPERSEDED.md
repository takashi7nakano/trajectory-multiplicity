# TRAJECTORY_MULTIPLICITY_PROTOCOL v1.0 — SUPERSEDED 記録

**記録日:** 2026-09-21
**対象:** `PROTOCOL_v1.0.md`（FROZEN 2026-09-21、本フォルダに未改変で保存）
**状態:** **SUPERSEDED BEFORE LESION OUTCOME COMPUTATION**
**根拠:** B REVIEW（2026-09-21）"PIPELINE HALT CONFIRMED. DO NOT RUN LESIONS."
**確定:** B VERDICT（2026-09-21）"ANALYTIC ASSESSMENT PASS. v1.0 SUPERSEDED."
**保存方針:** 永久保存。改変しない。

## 1. v1.0 は改変していない

`PROTOCOL_v1.0.md` は凍結時のまま保存する。本記録は別ファイルであり、v1.0 本文には
一文字も手を入れていない。v1.0 の凍結規則（§11 冒頭の freeze rule）はそのまま有効である。

## 2. supersede は outcome 計算の前に起きた

**registered removal seed は一つも instantiate していない。**
`data/tm_lesions_{member}_{ell}.npz` は一つも生成していない。
族の η_½ は一つも計算していない。

v1.0 が退けられた根拠は、すべて (i) 無摂動カーネルの汎関数、(ii) 族外トイ系での恒等式検証、
(iii) 解析的導出のいずれかである。**結果を見て endpoint を変えたのではない。**

これは §12 が記録する v0.1 → v0.2 → v0.3 の改訂と同じ性格のもので、
scientific null result ではなく protocol failure である。

## 3. supersede の理由（3 件、いずれも outcome 前に判明）

### 3.1 §5 P4 が族内定数（照会状 `2026-09-21_A_to_B_tm_P4_degenerate.md`）

層状ブロックグラフ上で §5 の指定どおり最大流を解くと、27 member × ℓ∈{4,8,16} すべてで
P4 = 1.000000000000。Σ_b μ^R_t(b) = 1 から各層のノード容量カットが値 1 に固定され、
反応流束自身がそれを達成するため、カーネルに依存しない。副次形 P4-int も
Birkhoff–von Neumann により n/ℓ で定数。

→ E1b・E2 が登録された形で計算不能。B 承認済み。

### 3.2 §5 P7 が全族で厳密計算不能（同上）

ブロック過程が強 lumpable でない member（Loc(σ) の ℓ≥2）では、ブロック過程は Markov 連鎖の
非 lumpable な関数であり系列エントロピーに閉形式がない。135 セル中 111 は厳密、
24 は Monte Carlo、うち 14 は 10⁵ 標本で推定量が飽和して値が定まらない。

→ P7 は現行の登録形では primary scalar たりえない。B 承認済み。

### 3.3 §4 primary configuration で endpoint 自体が退化（評価書 `2026-09-21_A_to_B_tm_analytic_assessment.md`）

persistent redirect の正確形は

    R(η)/R₀ = (1/n) Σ_{x∈S} P_x[T 歩のあいだ死行に到達しない]

であり、死行集合 Z(D) = ∅ のとき行和の正規化が committor の層内平坦性を保存するので
R(η)/R₀ = |S|/n となりカーネルに依存しない。死行は supp(P^R(x,·)) ⊆ D のときのみ生じるので、
**P^R(x,x) > 0 の member では原理的に生じない。**

族の対角最小値を確認した結果、**27 member 中 20 が狭義正の対角**を持つ。
それらでは全 η・全 ℓ で **η_½ = 1/2 が厳密**。変動しうるのは identity 以外の置換 7 member のみで、
すべて out-degree k = 1 であるため §3 C3 の ⟨k⟩ 整合部分族解析も実行できない。

→ 登録された介入は 𝒬₁ と 𝒬₂ を区別できない。カーネル依存性が support 性質のみを経由し、
二時刻同時分布が endpoint に入らないため。

## 4. 引き継がれるもの

| 項目 | 状態 |
|---|---|
| §2 模型クラス・埋め込み | 有効。T1–T5 で 27/27 通過（`tm_checks_v11.json`） |
| §3 族の生成規則 | 有効（member 数の記述誤り ≥30→27 を次版で訂正） |
| §7 seed 体系 | 有効。未消費 |
| §4 lesion family ℋ_ℓ の定義 | 有効。primary mode の選択は要再検討 |
| §5 P1–P3（𝒬₁） | 有効 |
| §5 P4, P7 | 撤回 |
| §5 P5, P6, P8 | exploratory として保持。**primary に昇格しない**（B 指示） |
| §6 E1a, E3, E4, E6 | 形式は有効。primary mode 変更の影響を要確認 |
| §6 E1b, E2, E5 | 再設計が必要 |
| §8 kill criteria | K1, K2, K5 有効。K3, K4, K6 は E1b/E2 再設計に依存 |

## 5. 次版の起草条件

B REVIEW の指示により、v1.1 の起草は本評価書の B レビュー後に行う。
A は v1.1 を起草していない。評価書 §5 に勧告を列挙したが、採否は B が決める。

v1.1 には以下を明記する。

- すべての改訂が lesion データを見る前に行われたこと
- 改訂の根拠が経験的性能ではなく導出であること
- §3 の member 数（≥30 → 27）
- 本記録へのポインタ

**更新（2026-09-21）:** v1.1 は凍結された（`PROTOCOL_v1.1.md`、SHA-256
4a34a4e1dea12adbb94bb5c067e1691a34d4301632f4600974ab890d2660a9dd、記録は
`PROTOCOL_v1.1_FREEZE_RECORD.md`）。B FINAL AUDIT PASS FOR FREEZE。凍結時点で outcome 未計算。

## 6. 関連文書

| ファイル | 内容 |
|---|---|
| `PROTOCOL_v1.0.md` | 凍結本文。未改変 |
| `2026-09-21_A_to_B_tm_P4_degenerate.md` | P4 退化の照会（B 承認済み） |
| `2026-09-21_A_to_B_tm_analytic_assessment.md` | 解析的プロトコル失敗評価 |
| `PHASE_1_2_REPORT_v11.md` | T1–T5 実残差、P1–P8 再実装 |
| `tm_checks_v11.json` | T1–T5 全残差 |
| `tm_predictors_v11.csv` | P1–P8 全値 |
| `tm_candidates_v11.csv` | 提案予測子 N_blk, S_blk |
| `_INVALID_haiku45_2026-09-21/` | 2026-09-21 Haiku 4.5 版の無効出力（隔離） |
