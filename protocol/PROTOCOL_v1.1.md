# TRAJECTORY_MULTIPLICITY_PROTOCOL v1.1 — FROZEN 2026-09-21

**Project:** Separate line from PHI_FLOW and from Paper I — *feed-forward reactive chain の
extensive-lesion robustness を再現するために、訪問ブロック数分布 P(V) の何次までの情報が必要か.*

**Status:** v1.1 **FROZEN 2026-09-21**。B FINAL AUDIT: **PASS FOR FREEZE**（2026-09-21）。
DRAFT → RC rev1（B REVIEW "MINOR REVISION"）→ RC rev2（B FINAL AUDIT "ONE BLOCKING
CLARIFICATION" への対応、member → family 集約規則の登録）を経て凍結。

**Freeze rule.** 本凍結以降、endpoint・閾値・族定義・予測子の役割・集約規則・kill criterion を
変更しない。凍結時点で outcome は何も計算されていない。**計算は §9 の順序に従い、PI の明示的な
指示があってはじめて開始する。**

**凍結宣言.** registered removal seed は未消費。`data/tm_lesions_*.npz`、`data/tm_curves_*.npz`、
`data/tm_PV_*.npz` はいずれも不存在。族の η_½・robustness curve・距離 D_{k,j} は一つも計算していない。
すべての改訂（v1.0 supersede を含む）は lesion outcome を見る前に行われ、根拠は経験的性能ではなく
(★) の恒等式・§4.2 の導出・登録質問の意味論である。集約規則を含め、curve 値を見て選んだ
登録要素は一つもない（§11.3、B FINAL AUDIT 確認済み）。
v1.0（FROZEN 2026-09-21）を supersede する。v1.0 本文は未改変で永久保存
（`PROTOCOL_v1.0.md`、記録は `PROTOCOL_v1.0_SUPERSEDED.md`）。

**Supersede の性格:** v1.0 は **SUPERSEDED BEFORE LESION OUTCOME COMPUTATION**。
registered removal seed は一つも instantiate されておらず、族の η_½ は一つも計算されていない。
v1.0 を退けた根拠はすべて (i) 無摂動カーネルの汎関数、(ii) 族外トイ系での恒等式検証、
(iii) 解析的導出である。**結果を見て endpoint を変えたのではない。**
scientific null result ではなく protocol failure である。

**経緯:** v1.0 §5 P4 の恒等的退化 → 照会 `2026-09-21_A_to_B_tm_P4_degenerate.md`（B 承認）
→ 解析評価 `2026-09-21_A_to_B_tm_analytic_assessment.md`（B PASS）
→ `REPAIR_SPEC_v1.1.md`（B MINOR REVISION、3 点修正のうえ本文起草へ）。§12。

**Scope:** feed-forward（層状、再訪なし）A→B Markov dynamics のみ。再帰的力学、連続状態系、
弱雑音極限、生物系、PHI_FLOW について何も主張しない。

**Executes:** Pipeline A。Pipeline B は凍結された族仕様からの独立再実装（A の結果凍結後）。

---

## 1. 登録質問と仮説

### 1.1 v1.0 からの変更の根拠

v1.0 §1 は「η_½ を担うのは 𝒬₁ か、𝒬₂ のうち構造量かスカラー量か」を問うた。
解析評価 §3.2 が、dead-end 介入について**除去集合を解析的に積分した正確な恒等式**

    E_D[R_dead(η)] / R₀ = E_ω[ C(n_b − V(ω), r) / C(n_b, r) ],   r = ⌈η n_b⌉     …… (★)

を与えた。V(ω) は軌道 ω が訪れる相異なるブロックの個数である。族外トイ系で独立検証済み（§7 T7）。

(★) は **robustness curve 全体が P(V) のみで決まる**ことを意味する。すなわち connectivity 情報は
この模型クラス・この介入において P(V) に射影される。したがって「connectivity か multiplicity か」は
分離可能な二分法ではなく、**識別可能な問いは情報階層のどの段で十分かに変わる**。

これは撤退ではない。outcome を見る前に、解析によって問いが厳密化された。

### 1.2 登録質問

> **Q.** dead-end robustness curve は、訪問ブロック数の**平均 E[V(ℓ)] だけ**で再現されるか。
> それとも訪問範囲分布 P(V(ℓ)) の**高次情報**を要するか。

情報階層:

    E[V]  ⊂  { E[V], Var(V) }  ⊂  P(V)

### 1.3 仮説

    H₀      登録解像度において E[V(ℓ)] のみで curve が再現される（closure-1 で十分）。
    H_dist  E[V(ℓ)] のみでは不十分。高次モーメントないし完全な P(V) が再現可能な情報を担う。

計算前の予想（A、2026-09-21）: **H_dist を予想する。**族内で CV = √Var(V)/E[V] が
0（並進対称な置換）から 1 超（Block(4,0.02) など）まで広がるため。
予想は記録であって endpoint ではない。

### 1.4 主張しない範囲

v1.0 §10 を継承する。加えて、(★) は **feed-forward・persistent block lesion・dead-end** の
下での恒等式であり、他の介入・他の力学クラスに一般化しない。

---

## 2. 模型クラスと構成規則（v1.0 §2 から変更なし）

**鎖.** 層 t = 0..T、各層 n 位置、x ∈ {0..n−1}。層 T の後に吸収 B、全層から到達可能な吸収 F。
A = 層 0、μ_A は位置上一様。feed-forward、再訪なし。

**構成規則（B11）.** 反応核を先に指定してから埋め込む。

1. 反応周辺分布 μ^R_t を選ぶ。primary family: 一様、μ^R_t(x) = 1/n。
2. 反応核 P^R_t を §3 のメニューから選ぶ（二重確率）。
3. 一様 leak で埋め込む: P_t(x,y) = (1−λ) P^R_t(x,y)、P_t(x,F) = λ。
   このとき q_t(x) = (1−λ)^{T−t}、R₀ = (1−λ)^T、Doob 条件付き過程は厳密に P^R_t。
   secondary embedding（報告のみ）: 位置依存 leak λ_t(x)。

**サイズ.** n = 256、T = 64、λ = 1 − 0.5^{1/64} = 0.010771986806、R₀ = 0.500000000000。
解像度 ℓ ∈ {1, 2, 4, 8, 16}、ブロック index b = ⌊x/ℓ⌋、n_b = n/ℓ。
分割原点感度チェック（C1）は ℓ ≥ 4 で s ∈ {ℓ/4, ℓ/2, 3ℓ/4}。

T1–T5 は v1.0 の実装で 27/27 member 通過済み（`tm_checks_v11.json`、最大残差 2.1e-12）。

---

## 3. 族と条件

### 3.1 メニュー（v1.0 §3 から変更なし。member 数の記述誤りのみ訂正）

すべて n = 256 上の二重確率核、周期境界、時間一様。

| 種別 | 個数 | パラメータ |
|---|---|---|
| Loc(σ) | 6 | σ ∈ {1, 2, 4, 8, 16, 32} |
| Perm(π) | 8 | identity、cyclic shift s ∈ {1, 4, 16, 64}、固定乱数置換 3（seed 2026101–2026103） |
| Mix(c) | 5 | c ∈ {0.02, 0.05, 0.1, 0.2, 0.5} |
| Block(M,c) | 6 | M ∈ {4, 16, 64} × c ∈ {0.02, 0.1} |
| ShiftMix(c) | 2 | s = 16、c ∈ {0.05, 0.2} |
| **合計** | **27** | |

**訂正:** v1.0 §3 は「the menus below give ≥ 30」と記したが、メニューの実数は **27** である。
登録目標 N_fam ≥ 24 は満たす。

核は厳密構成とし、Sinkhorn 反復を用いない（既に二重確率な行列への反復は浮動小数点雑音を
注入するだけである）。T1 の行和・列和残差は 0〜1.1e-14。

### 3.2 条件

| | 内容 | v1.0 → v1.1 |
|---|---|---|
| **C1** | 全予測子がブロックラベル置換 10 種で不変（1e-10） | **維持** |
| **C2** | 一様反応周辺分布（逸脱 > 1e-10 の member は除外、K2） | **維持**。dead-end では endpoint を潰さない |
| **C3** | 次数・流束整合部分族 | **置換**（下記 3.3） |
| **C4** | 対照 K-a, K-b | **維持・強化**（下記 3.4） |
| **C5** | カット構造整合（Loc vs Mix の min-cut 5% 一致） | **撤回**。min-cut は退化した P4 の量 |
| **C6** | 同一 coherence・異なる配線 | **維持。strong control であって primary endpoint ではない**（下記 3.5） |

### 3.3 C3 の再設計 — 事後の matched subfamily を作らない

v1.0 C3 は ⟨k⟩ 整合部分族での再解析を求めたが、(★) の下で robustness を決めるのは ⟨k⟩ ではなく
V の法である。また E[V] は既に全族で計算済みであるため、**事後に「よく合う対」を選ぶ自由度を
持ち込まない。**

**C3′:** 全 27 member を用い、**E[V(ℓ)] を連続座標として扱う。**別途の matching は行わない。

理由: **closure 階層そのものが条件付き解析になっている。**closure-2 は closure-1 と
**構成上同じ平均を持つ**（§5.3）ので、「Var(V) が E[V] を超える情報を足すか」は
E[V] を固定したうえでの比較である。人工的な部分族は不要である。

### 3.4 C4 対照（データ取得前に実施）

| | 内容 | v1.1 での位置づけ |
|---|---|---|
| **K-a** | Loc(1) に独立補助座標 z ∈ {0..15} を毎層一様再抽出で付加。lesion は x のみに定義 | V は x のブロックで定義されるので z は V を変えない。**E[V], Var(V) が厳密に一致すること（1e-12）を検証する** |
| **K-b** | Loc(1) の層 T/2 に二峰待機ゲート（0 または 8 層、等確率）。恒等核の待機層 8 枚を確率 1/2 で通過 | 恒等待機層は新しいブロックを訪れないので V は不変。**E[V], Var(V) が gate なしの Loc(1) と厳密に一致すること（1e-12）を検証する** |

**v1.0 からの改善:** v1.0 では両対照を P7（itinerary N_eff）で検証する設計だったが、
P7 は Loc 系で Monte Carlo 飽和するため厳密検証が不可能だった。V を用いる v1.1 では
**両対照とも厳密に検証できる。**

raw path entropy は K-a で T ln 16 だけずれるため使用しない（v1.0 と同じ）。

### 3.5 C6 — strong control

対 {identity, shift(16), random permutation} は ℓ=1 で C = T ln n = 354.891356 が共通、
対 {Mix(c), ShiftMix(c)} は同 c で C が共通（いずれも A 実装で厳密に再現済み）。
一方 E[V] は大きく異なる（ℓ=16 で identity 1.0000、shift-16 16.0000、random-0 14.4297）。

**C6 は「旧スカラー coherence C が新しい十分座標の階層を決定しない」ことの
sanity / strong control である。primary endpoint には昇格させない。**
primary question は closure-1 と truth の比較である（§6.2）。

---

## 4. Lesion family と介入

### 4.1 Footprint と除去集合（v1.0 §4 から変更なし）

ℋ_ℓ の元はブロック index b ∈ {0..n_b−1}。強度 η の extensive lesion は
r(η) = ⌈η n_b⌉ 個のブロック index の一様無作為部分集合。
時間構造 — **primary: persistent**（除去ブロックは全層で除去）。secondary: 層ごと独立除去。

### 4.2 Mode — primary を dead-end に反転（D3 の反転）

| | v1.0 | v1.1 |
|---|---|---|
| primary | redirect | **dead-end** |
| secondary | dead-end | redirect |

根拠（解析評価 §1.3, §3.1、B 承認済み）:

    redirect:  R(D)/R₀ = (1/n) Σ_{x∈S} P_x[T 歩のあいだ死行に到達しない]
    dead-end:  R(D)/R₀ = P[x₀,…,x_T がすべて S に留まる]

redirect では死行集合 Z(D) = ∅ のとき行和の正規化が committor の層内平坦性を保存するので
**R(D)/R₀ = |S|/n となりカーネルに依存しない。**死行は supp(P^R(x,·)) ⊆ D のときのみ生じ、
**P^R(x,x) > 0 の member では原理的に生じない。**族の 27 member 中 **20 が狭義正の対角**を持ち、
それらでは全 η・全 ℓ で **η_½ = 1/2 が厳密**である。変動しうるのは identity 以外の置換 7 member のみで、
すべて out-degree k = 1 である。

dead-end は全 member で非退化であり、(★) によって curve 全体が P(V) で決まる。

**redirect の新しい位置づけ.** v1.0 E5（redirection gain）は実験的 endpoint だったが、gain は解析的に既知:

    G(η)/R₀ = |S|/n − G_V(1−η)   （対角が正の 20 member）
    G(η)/R₀ = 0                   （identity 以外の置換 7 member、redirect ≡ dead-end）

gain が正になるのは endpoint が定数になる member に限られるため、
**E5 は endpoint から外し、恒等式検証項目（§7 T7）に移す。**

### 4.3 η グリッド — 解析的極限のみから導出

**特定 member の計算結果は一切用いていない。**

(★) より R/R₀ = E[(1−η)^V]、V ∈ [1, n_b] なので η_½ = 1 − 2^{−1/V}。解析的 bracket:

| 極限 | η_½ |
|---|---|
| V = 1（完全閉じ込め） | **0.500000** |
| V = n_b（最大拡散） ℓ=16 | 0.042397 |
| 〃 ℓ=8 | 0.021428 |
| 〃 ℓ=4 | 0.010772 |
| 〃 ℓ=2 | 0.005401 |
| 〃 ℓ=1 | **0.002704** |

v1.0 のグリッド {0.05,…,0.95} は全解像度で下側極限を bracket しない。

**登録グリッド.** u = −ln(1−η) について幾何的に刻む（η_½ ↔ u = ln2/V）。

    u_k = ln 2 · 2^{−k/4},   k = 0, 1, …, 4·log₂ n_b,   η_k = 1 − exp(−u_k)

1 オクターブあたり 4 点。ℓ=16 で 17 点、ℓ=1 で 33 点。両端が上表の解析的極限に一致する。

    ℓ=16: 0.5000, 0.4417, 0.3875, 0.3378, 0.2929, …, 0.05943, 0.05021, 0.04240
    ℓ=1 : 0.5000, 0.4417, 0.3875, 0.3378, 0.2929, …, 0.00382, 0.00321, 0.00270

**補間誤差の事前上限.** R(u) = E[e^{−uV}] は完全単調、R'' = E[V²e^{−uV}] > 0。
刻み比 2^{−1/4} で Δu/u = 0.1591。半値点 u ≈ ln2/E[V] における線形補間誤差は

    |ε_R| ≤ (1/8)(Δu)² E[V²] = 1.52×10⁻³ · (1 + CV²)

**登録上限 |ε_R| ≤ 5×10⁻³。**CV² > 2.3 の member は当該 ℓ で 1 オクターブあたり 8 点に倍化する
（規則として事前登録する。member ごとの判定は既計算の CV から決まり、outcome を見ない）。

> **重要（B 指示）:** この 5×10⁻³ は **R の線形補間誤差の上限**である。
> **η_½ に対する誤差保証ではない。**η_½ への読み替えを行わない。
> §6 の判定統計量は R の曲線距離であってグリッド点上で評価するため、補間を経由しない。

---

## 5. 予測子と推定量

### 5.1 撤回・保持

| | v1.0 | v1.1 | 理由 |
|---|---|---|---|
| P4（fractional max-flow）| 𝒬₂ structural primary | **撤回** | 27 member × ℓ∈{4,8,16} すべてで 1.000000000000。Σ_b μ^R(b)=1 から各層のノード容量カットが 1 に固定され、反応流束自身がそれを達成する。副次形 P4-int も Birkhoff–von Neumann により n/ℓ で定数 |
| P7（itinerary N_eff）| 𝒬₂ scalar primary | **撤回** | 非 lumpable member（Loc(σ), ℓ≥2）で系列エントロピーに閉形式がない。135 セル中 111 厳密、24 が MC、うち 14 が 10⁵ 標本で飽和 |
| P1–P3 | 𝒬₁ | **保持**（報告のみ。一様族では構成上定数） | |
| P5, P6, P8 | exploratory | **保持。primary に昇格しない**（B 指示） | |

**新しい structural graph statistic は探さない。**(★) により endpoint は P(V) のみで決まる。

### 5.2 登録予測子

| 記号 | 定義 | 役割 | 計算 |
|---|---|---|---|
| **P_scalar-1** | N_blk(ℓ) = E[V(ℓ)] | 一次。H₀ の内容 | **厳密**。E[V] = Σ_b (1 − a_b)、a_b = P[ブロック b を一度も訪れない] = (1/n)·𝟙ᵀ(P^R\|_{bᶜ})^T 𝟙。n_b 回の T ステップ再帰 |
| **P_scalar-2** | S_blk(ℓ) = Var(V(ℓ)) | 二次。registered second-order predictor | **厳密**。単一・対ブロック回避から。n_b(n_b−1)/2 回の追加再帰 |
| **P_oracle** | 完全な P(V(ℓ)) | **解析的 ground truth benchmark。競争する predictor ではない** | §5.5 |

Var(V) は v1.0 相当の「structural predictor」ではない。**同一の一次元確率変数 V のスカラー要約**であり、
pairwise co-visitation が計算に現れることは variance を connectivity statistic にしない（B 指摘、A 受諾）。

計算済み（`tm_candidates_v11.csv`、SHA-256 268b95ab1635fb452cd5636c081b250e0429fa6d35d7cf6d395c632f400e61c1）。
ℓ=4 で N_blk は 1.0–64.0（27 member 全て相異なる値）、S_blk は 0–126.0。

### 5.3 Closure — 対称な最大エントロピー構成

support {1, …, n_b} 上で、同一の最大エントロピー原理により定義する。

    **closure-1   p₁(v) ∝ exp(λ₁ v)              平均のみを拘束**
    **closure-2   p₂(v) ∝ exp(λ₁ v + λ₂ v²)      平均と分散を拘束**

Lagrange 乗数は登録モーメントに一意に合わせる。両段の違いは**利用可能な情報量だけ**であり、
closure 構成の違いが交絡しない。

**退化の登録規則.**

- closure-1: E[V] = 1 または n_b のとき端点の点質量（Perm(identity) が該当）。
- closure-2: Var(V) = 0 のとき E[V] における点質量。

**well-posedness 確認済み**（既計算の無摂動モーメントのみ使用、lesion なし）。全 5 解像度で
22 member が通常解、5 member が Var(V)=0 につき退化（並進対称な置換: identity, shift-1, shift-4,
shift-16, shift-64）、モーメント空間外・不収束は 0。退化 25 セルの E[V] は**すべて厳密に整数**
（残差 0.00e+00）であり、点質量の位置に丸めの任意性がない。

### 5.4 Primary estimator — 除去集合を解析的に積分する

(★) の右辺は**経路 ω についてのみ**の期待値であり、除去集合 D は閉形式で積分済みである。したがって

    **R̂(η)/R₀ = (1/M) Σ_{i=1}^{M} w(V_i, r(η)),   w(V, r) = C(n_b − V, r) / C(n_b, r)**

（n_b − V < r のとき w = 0）。V_i は無摂動反応過程から抽出した M 本の軌道の訪問ブロック数。

性質:

- **removal-set noise が構成上ゼロ。**標本数で小さくするのではなく解析的に除去される。
- 一組の軌道標本から**全 η グリッドの curve が同時に**得られる。
- **一組の軌道標本から全解像度の V が同時に**得られる（軌道は ℓ に依らない）。解像度間が完全に結合する。
- lumpability も閉形式も不要な**単一の汎用算法**であり、Pipeline B が独立再実装しやすい。
- 残る誤差は [0,1] 有界汎関数の軌道標本誤差のみ。supp(V) ⊆ {1,…,n_b} と小さい。

**登録標本数 M = 10⁶。**推定量は [0,1] 有界変数の平均なので、グリッド点ごとに SE ≤ 1/(2√M) = 5.0×10⁻⁴。

**登録軌道 seed:** 2026400000 + 10⁴ · (member index)。解像度に依存しない（同一標本を全 ℓ で使う）。
v1.0 の removal-set seed 体系（2026200000 + 10⁴·member + 10²·ℓ + η index）は**変更しない**
（B 指示）。§7 の validation で用いる。

### 5.5 Oracle — 厳密 P(V) が得られる 13 member

| 経路 | 対象 | 算法 |
|---|---|---|
| **A 直接列挙** | 置換 8 member | 決定論的核。n 個の始点それぞれの軌道を辿る。O(nT) |
| **B 占有数閉形式** | Mix(c) 5 member | 跳躍 K〜Binom(T,c)、着地点 iid 一様。m = K+1 回の一様配置の占有数分布。O(T·n_b²) |

    P(V=v) = Σ_k Binom(T,k;c) · C(n_b,v) · Σ_{j=0}^{v} (−1)^j C(v,j) ((v−j)/n_b)^{m},  m = k+1

両経路とも独立に計算した厳密 E[V]（ブロック回避再帰）と照合済み（差 0.0e+00 〜 2.9e-10、§7 T8）。

**oracle は scientific endpoint ではなく implementation validation である。**
**oracle を 13 member より広げるための追加の厳密算法は実装しない**（B 指示）。
Loc のレンジ DP、Block / ShiftMix の厳密算法探索は行わない。

**厳密 P(V) の全族計算は実行不可能である。**訪問ブロック集合の格子上の包除が必要で 2^{n_b} 部分集合を要する
（ℓ=1 で 2²⁵⁶）。一般の Markov 連鎖について訪問セル数分布を多項式時間で与える算法は知られていない。

### 5.6 P7 が破綻し V が成立する理由

**supp(V) ⊆ {1, …, n_b}。**P7 の itinerary 分布は台が ~n_b^T 個の原子を持ち、10⁵ 標本で
135 セル中 14 が飽和した。P(V) は高々 n_b 個のビンのヒストグラムであり標本抽出で正確に推定できる。
同じ手法が P7 で破綻し P(V) で成立する理由は**台の大きさ**である。

実測（ℓ=8、n_b=32、2×10⁵ 標本）で標本 Ê[V] は厳密 E[V] と全 member で SE 内一致（§7 T8）。

---

## 6. Endpoint と判定規則

### 6.1 曲線距離（member ごとに定義される判定統計量）

判定量は **member ごと**に定義される。member を j = 1..27、その無摂動核を P^R_j とする。
member j の closure-k curve は、member j 自身の（厳密な）モーメントから作った p_{k,j} を
(★) の重み w(v, r) = C(n_b − v, r)/C(n_b, r) に代入して定義する。

    R_{k,j}(η)/R₀ = Σ_{v=1}^{n_b} p_{k,j}(v) · w(v, r(η)),   k = 1, 2

member j の truth R*_j(η) は、oracle 13 member では厳密（§5.5）、残り 14 member では §5.4 の
D 積分済み推定量。member ごとの曲線距離:

    **D_{k,j}(ℓ) = max over registered grid η of | R_{k,j}(η) − R*_j(η) | / R₀**

グリッド点上で評価するので**補間を経由しない**（§4.3 の注記）。

**誤差予算（member ごと）:** oracle member（13）は誤差 0。推定 member（14）はグリッド点ごと
SE ≤ 5.0×10⁻⁴、グリッド点数 ≤ 33 の sup に対する Bonferroni 調整込みで ≤ 2.0×10⁻³。
closure curve R_{k,j} は**全 27 member で厳密**である（厳密なモーメント E[V]_j, Var(V)_j から構成するため）。
したがって D_{k,j} の推定誤差は R*_j の側のみに由来し、oracle 13 member では 0 である。

### 6.2 E2 — primary endpoint（member → family の集約規則を明示）

**登録閾値 τ_R = 1.0×10⁻²**（誤差予算 2.0×10⁻³ の 5 倍）。

#### 6.2.1 member ごとの判定

member j について、登録解像度 ℓ ∈ {1, 2, 4, 8, 16} のうち closure-1 が許容を超える個数を

    **m_j = #{ ℓ ∈ {1,2,4,8,16} : D_{1,j}(ℓ) > τ_R }**

と置く。member ごとの判定:

    **m_j = 0**   → member j は平均で再現可能（H₀-consistent）
    **m_j ≥ 2**   → member j は平均を超える情報を要する（H_dist member）
    **m_j = 1**   → member j は inconclusive（解像度に限定）

三つは網羅的かつ排他的。**同一 member が 2 解像度以上で超過する**ことを H_dist member の
要件とする（解像度をまたいだ pool ではない）。

#### 6.2.2 family への集約規則（existential rule。curve を見る前に登録する）

family verdict は 27 member の {m_j} の**決定論的関数**として次のとおり定める。

    **family H₀**          ⟺ すべての member で m_j = 0
    **family H_dist**      ⟺ 少なくとも一つの member で m_j ≥ 2
    **family inconclusive** ⟺ 上のいずれでもない（m_j ≥ 2 の member はないが、m_j = 1 の member がある）

三帰結は {m_j} 上で網羅的かつ排他的。**τ_R 以外に新たな自由パラメータを導入しない**
（member 集約の重み・割合・順位・member 選択のいずれも用いない）。

**規則の根拠（研究質問との対応）.** 登録した座標充足性の主張——情報階層
E[V] ⊂ {E[V], Var(V)} ⊂ P(V)——は「E[V] が robustness の**十分統計量**か」を問う。
充足性は普遍的性質である：E[V] が family にとって十分であるとは、**すべての member の
curve を各自の平均だけで再現できる**ことに他ならない。したがって、ある一つの member の
curve が（2 解像度以上で頑健に）平均を超える情報を要するなら、その一点で
「E[V] は family にとって十分」は反証される。これが existential rule を選ぶ理由である。
member をまたいだ多数決・平均・中央値・割合はいずれも「典型的 member では十分か」という
**別の質問**（B が「科学的意味が異なる」と指摘したもの）に答えるものであり、
情報階層の枠組みが問うているのは充足性＝普遍性である。

#### 6.2.3 H_dist 内の分類（curve を見た後）

D₂ は H_dist の**必要条件ではない**（B REVIEW 2026-09-21）。H_dist member それぞれについて、
D_{1,j}(ℓ) > τ_R となった各 ℓ で

      **D_{2,j}(ℓ) < D_{1,j}(ℓ)/2**  → その ℓ で **H_dist-2nd**：二次モーメントが近似を実質的に改善
      **それ以外**                    → その ℓ で **H_dist-higher**：**二次モーメントを超える情報を要する**

**H_dist-higher は H_dist の失敗ではない。**情報階層の第三段（P(V) のより高次の情報が必要）を
指す所見として報告する。oracle 13 member では完全な P(V) との差 D_P が定義上 0 であることを併記する。

#### 6.2.4 主要報告対象

**primary reported object は 27 member × 5 解像度の表**：各セルの D_{1,j}(ℓ), D_{2,j}(ℓ)、
各 member の m_j と分類。family H₀ / H_dist / inconclusive の verdict はその**導出要約**である。

**H₀ の登録文言**（family H₀ のときのみ）:
「登録された族・登録されたサイズ・登録された解像度において、族の**全 member**について
dead-end robustness curve が訪問ブロック数の平均のみで登録許容内で再現された」。一般の主張ではない。

**inconclusive の登録文言**（family inconclusive のとき）:
「平均を 2 解像度以上で超えた member は無く、単一解像度でのみ超えた member があった。
登録された多解像度基準を満たさないため、仮説を選択せず、解像度に限定された所見として報告する」。
**H₀ 側として報告しない。**

#### 6.2.5 推定雑音に対する安全性（M から今言える）

closure curve は全 27 member で厳密。truth は oracle 13 member で厳密、他 14 member で
グリッド点ごと SE ≤ 5.0×10⁻⁴。τ_R = 1.0×10⁻² は sup 誤差予算（≤ 2.0×10⁻³）の 5 倍であり、
真に平均再現可能な member の D̂_{1,j} が τ_R を一解像度でも超える確率は、点あたり ~20σ の事象で
無視できる（member 数 ×グリッド点数を掛けても同様）。**推定雑音は H_dist を捏造できない。**
さらに、V が集中する member（置換 8、Mix 5）は**すべて oracle 集合**にあり誤差 0 であるため、
family verdict の H_dist driver は厳密 curve に基づき、推定に依存しない。

#### 6.2.6 規則の帰結の透明性（モーメントのみ。curve は計算しない）

closure-1 は support {1,…,n_b} 上の平均のみ最大エントロピー分布であり、E[V] が support の
内部にある限り**退化しない（広がる）**。点質量になるのは E[V] ∈ {1, n_b} のときに限る（§5.3）。
族は構成上（§3）、訪問ブロック数が内部の平均に**集中**する member（Var(V) = 0、既計算）を含む。
それらでは closure-1 が集中した truth を再現できないため、**family verdict は H_dist が予期される**。
これは既計算のモーメントから導かれる**定性的**帰結であり、透明性のために開示する。
**距離 D_{k,j} 自体は凍結前に計算しない。**

真に未計算で開いているのは、広がった分布を持つ member（Loc, Block, ShiftMix）の
D_{1,j} / D_{2,j} の**大きさ**と、族全体での **H_dist-2nd / H_dist-higher の分類**である——
すなわち、どの member・どの regime で二次モーメントが十分か、二次を超える情報を要する member が
存在するか。これが本 endpoint の科学的中身であり、curve を見るまで決まらない。

### 6.3 その他の endpoint

| | 内容 | v1.0 → v1.1 |
|---|---|---|
| **E1** | 族内変動の記述的診断（下記 6.4） | **降格**（B 指示）。**decision weight なし。閾値も判定もない。いかなる kill criterion にも入らない** |
| **E2** | §6.2 | **新 primary** |
| **E1b** | Δρ₁ = \|ρ(P4)\|−\|ρ(P2)\| | **撤回**（P4 撤回） |
| **E3** | C6 対の curve 差 | **strong control として維持。primary へ昇格しない**（B 指示）。旧スカラー coherence C が新しい十分座標の階層を決定しないことの sanity |
| **E4** | C5 残差 | **撤回**（C5 撤回） |
| **E5** | redirection gain | **撤回**。§4.2 のとおり解析的に既知。§7 T7 の恒等式検証へ |
| **E6** | 層ごと独立除去での時間構造チェック | **維持**。全 member で生存が (1−η)^T になることの確認 |
| **E7** | 探索的 out-of-sample | **縮小して維持。E[V] / Var(V) / P(V) の階層内に限定する。新しい統計量は探さない**（B 指示）。分割 seed 2026300 |

### 6.4 E1 — 族内変動の記述的診断（decision weight なし）

v1.0 の E1a（excess variance + 階層 bootstrap）は、各 member の η̂_½ が removal-set 雑音を
持つことに由来する構成だった。§5.4 で除去集合を解析的に積分したためその雑音は存在しない。
残るのは軌道標本誤差のみであり、oracle 13 member ではそれも 0 である。

B REVIEW（2026-09-21）の指示により、**E1 は記述的診断に降格する。**

報告する内容（閾値なし、判定なし）:

- 各グリッド点 η における R*(η)/R₀ の member 間の範囲・四分位・全 27 member の値
- 同じグリッド点における member ごとの推定誤差（oracle 13 member は 0、他は SE ≤ 5.0×10⁻⁴）
- 記述的要約としての η̂_½ の member 間分布（§4.3 の注記どおり、判定には用いない）

**E1 はいかなる仮説の宣言にも、いかなる kill criterion にも寄与しない。**
族に変動があるかどうかは §6.2 の判定に影響しない。

報告規則: 全 member、全 ℓ、両 embedding、両 mode、両時間構造を、endpoint の結果によらず報告する。

---

## 7. Code checks、seed、凍結出力

### 7.1 検査

| | 内容 | 許容 | 状態 |
|---|---|---|---|
| **T1** | 二重確率性と μ_A の T 層伝播 | 1e-12 | **27/27 済**（最大 1.1e-14） |
| **T2** | 埋め込み鎖の後退再帰 committor vs (1−λ)^{T−t}、Doob 変換 vs P^R、R₀ = 1/2 | 1e-12 | **27/27 済**（最大 1.4e-14） |
| **T3** | H(Γ\|𝒮) = Σ_t H_t − C（3 経路を独立構成、全 ℓ） | 1e-10 | **27/27 済**（最大 2.1e-12） |
| **T4** | single-edge lesion 100/member、解析式 vs 総当たり後退再帰 | 1e-12 | **27/27 済**（最大 1.3e-14） |
| **T5** | C1 ブロックラベル置換 10 種での全予測子不変性 | 1e-10 | **27/27 済**（最大 2.8e-13） |
| **T6** | 対照 K-a, K-b（§3.4）。E[V], Var(V) が厳密一致 | 1e-12 | **未実施** |
| **T7** | 解析恒等式 4 本を族外トイ系（n=8, T=4）で総当たり検証: redirect 作用素形、dead-end 作用素形、Z(D)=∅ での \|S\|/n、(★) | 1e-12 | **済**（最大 1.2e-14） |
| **T8** | oracle 13 member で厳密 P(V) と §5.4 推定量が一致 | **グリッド点ごと 4 SE**（下記注） | **部分実施**（E[V] 照合済み。curve 照合は未実施） |
| **T9** | closure-1, closure-2 が目標モーメントを再現 | 1e-8 / 1e-6 | **済**（27 member × 5 ℓ、失敗 0） |

**T8 の許容について（B 指示）:** 4 SE は**グリッド点ごとの実装許容**であり、
**曲線全体についての同時信頼記述ではない。**グリッド点数 ≤ 33 に対する多重性の調整を含まず、
統計的推論ではなく実装の一致確認として用いる。

**K1 の拡張:** T8 の不一致は code failure として K1 の対象とする。

### 7.2 Seed（v1.0 体系を維持。B 指示により変更しない）

    族の乱数置換      2026101–2026103
    removal set       2026200000 + 10⁴·(member index) + 10²·(ℓ index) + (η index)   ※ validation 用
    E7 分割           2026300
    軌道標本（新規）   2026400000 + 10⁴·(member index)

**開示（B 裁定: contamination 扱いしない、開示は維持、seed 変更不要）.**
実行可能性確認のため、ℓ=8 の 7 member について P̂(V) のヒストグラムを一時的に生成した
（seed 4242、2×10⁵ 標本）。R(η)、η_½、removal set はいずれも形成せず、ヒストグラムは
ディスクに書き出していない。出力したのは E[V]、台の大きさ、SE のみ。
ただし P(V) は curve の十分統計量であるため、原理的にはこの 7 セルについて curve を決定しうる
情報を A は一時的に保持した。seed 4242 は登録系列（2026400000 系）に属さないため衝突はない。

### 7.3 凍結出力

    data/tm_family.json                  核仕様（27 member）
    data/tm_moments_{ell}.csv            E[V], Var(V)（凍結、SHA-256 を manifest に）
    data/tm_PV_{member}.npz              推定 P̂(V)（全 ℓ、同一軌道標本）
    data/tm_PV_exact_{member}.npz        厳密 P(V)（oracle 13 member）
    data/tm_curves_{member}_{ell}.npz    R*, R₁, R₂
    data/tm_lesions_{member}_{ell}.npz   removal-set validation
    results/tm_manifest_v11.md           全 SHA-256

---

## 8. Kill criteria（いずれか一つで null report、再解釈なし）

| | 内容 | v1.0 → v1.1 |
|---|---|---|
| **K1** | T1–T9 のいずれかが失敗 → pipeline invalid。修正し、版を上げ、族データを見る前に再検査 | **維持・拡張**（T7–T9 を追加。oracle 不一致を含む） |
| **K2** | C2 が破れる member（周辺分布逸脱 > 1e-10）→ 当該 member を除外 | **維持** |
| **K3** | **family H₀**（§6.2.2、全 member で m_j = 0）→ 登録族において全 member が平均で十分。§6.2.4 の限定文言で報告し、線を止める | **置換**（新 family H₀ に読み替え） |
| **K4** | **family inconclusive**（§6.2.2、m_j ≥ 2 の member なし、m_j = 1 の member あり）→ **仮説を選択しない。**§6.2.4 の限定文言で報告する | **書き換え**（B 指示）。v1.0／v1.1 DRAFT の「H₀ 側を報告」は誤り。**H₀ 側として報告しない** |
| **K5** | 予測子が C1（T5）に失敗 → 当該予測子を除外、研究は継続 | **維持** |
| **K6** | primary embedding と secondary embedding が §6.2.2 の family verdict のうち**異なるもの**を与える（H₀ / H_dist / inconclusive のいずれか 2 つ）→ embedding 依存と宣言し、仮説を選択しない | **書き換え**（E1b 撤回、および §6.2 の family 集約規則に伴う） |

K3, K4, K6 はいずれも **§6.2.2 の {m_j} → family verdict 集約規則**にのみ依拠し、
単一の家族レベル m には依拠しない（B FINAL AUDIT 2026-09-21 の指摘に対応）。
**family H_dist**（∃j: m_j ≥ 2）は期待される結果であって研究の失敗ではない（旧 K4 の趣旨を継承）。

**E1 は kill criterion に寄与しない**（§6.4、B 指示）。

---

## 9. 実行順序

    コードと seed の凍結
      → T1–T5（済）
      → T6（対照 K-a, K-b）、T7（済）、T9（済）
      → E[V], Var(V) を全 member・全 ℓ で凍結（manifest）
      → 軌道標本 M = 10⁶ を member ごとに抽出、全 ℓ の P̂(V) を凍結
      → oracle 13 member の厳密 P(V) を計算、T8 で推定量と照合
      → R*, R₁, R₂ を構成
      → E2（primary）→ E3（strong control）→ E6 → E7（探索）
      → E1（記述的診断。decision weight なし）
      → secondary（redirect mode、層ごと独立除去、secondary embedding）
      → removal-set validation（v1.0 seed 体系、N_rem = 200）
      → Pipeline A 結果を凍結
      → Pipeline B が data/tm_family.json から独立再実装。
        PASS 条件: E[V], Var(V) が 1e-8 で一致し、R* が推定誤差内で一致

**lesion / outcome computation は v1.1 が B レビューを受け明示的に凍結されるまで行わない。**

---

## 10. 主張しない範囲

v1.0 §10 を継承する。再帰的力学（循環 toy が先に必要）、連続状態系、カオス系、弱雑音系、
生物系・脳・同一性・意識、PHI_FLOW、Paper I について何も主張しない。
「並列経路が robustness を生む」は既知であり主張しない。

加えて v1.1 固有:

- (★) は feed-forward・persistent block lesion・dead-end の下での恒等式であり、
  他の介入・他の力学クラスに一般化しない。
- H₀ / H_dist の結論は**登録された族・サイズ・解像度**についてのものであり、
  「一時刻統計が extensive-lesion robustness を決める」といった一般命題ではない。
- v1.0 の 𝒬₁/𝒬₂ 階層は、(★) の射影が成り立つ範囲でのみ情報階層
  E[V] ⊂ {E[V], Var(V)} ⊂ P(V) に置き換わる。模型クラス外への持ち出しはしない。

---

## 11. DRAFT に対する B REVIEW の処理

B REVIEW OF PROTOCOL v1.1 DRAFT（2026-09-21）: **MINOR REVISION BEFORE FREEZE.**

### 11.1 承認された項目（変更なし）

| DRAFT §11 の照会 | B の裁定 |
|---|---|
| 1. M = 10⁶、軌道 seed 系列 2026400000 | **承認** |
| 3. 判定統計量を sup ノルム曲線距離とし、η_½ は記述的要約に留める | **承認** |
| 4. T8 の一致許容 4 SE | **承認**（ただし実装許容であって同時信頼記述ではないと明記すること → §7.1 に反映） |
| 6. C3′（全族 + E[V] 連続座標、matching なし） | **B の指示を満たしている** |
| 2 のうち τ_R = 1.0×10⁻² | **承認** |

### 11.2 修正した 4 点

| B の指示 | 反映先 |
|---|---|
| D₂ < D₁/2 は **H_dist の必要条件ではない**。closure-1 の失敗（D₁ > τ_R が 2 解像度以上）だけで H_dist が成立する。D₂ は H_dist 内の**分類**に用いる。条件を満たさない場合は「二次モーメントを超える情報が必要」と報告し、**H_dist の失敗として分類しない** | **§6.2** を三帰結（m=0 / m≥2 / m=1）に書き換え、H_dist-2nd / H_dist-higher の分類を追加 |
| K4 を改訂。H₀ でも登録された多解像度 H_dist 基準でもない場合は **inconclusive / resolution-limited** として報告し、**「H₀ 側」としない** | **§8 K4** を m = 1 の場合に書き換え。§6.2 に登録文言を追加 |
| E1 は仕様不足。**記述的な族内変動診断に降格**し decision weight を持たせない（推奨）。さもなくば凍結前に統計量と閾値を完全に指定する | **§6.3 / 新設 §6.4**。降格を採用。閾値なし、判定なし、kill criterion に寄与しない |
| T8 に「4 SE はグリッド点ごとの実装許容」と明記 | **§7.1** の許容欄と直後の注記 |

あわせて §8 K6 を、§6.2 の三帰結構造に合わせて書き換えた（embedding 間で**異なる帰結**を与える場合）。

（rev1 で §6.2 に置いた三帰結 m=0 / m≥2 / m=1 は**単一の家族レベル m** を前提にしていた。
rev2 でこれを member 解決型に置き換えた。次項参照。）

### 11.3 B FINAL AUDIT（2026-09-21）: member → family 集約規則の登録

B FINAL AUDIT: **ONE BLOCKING CLARIFICATION.**「ご照会の 2 点（rev1 の三帰結の網羅性、
E1 降格）は PASS。ただし §6.1 の距離が member ごとに D_{k,j}(ℓ) と定義される一方、family verdict
への集約規則が未登録であり、curve を見た後に選べない。E2 と K3/K4/K6 のみを更新して集約規則を
完全に指定せよ。」

**未登録だった自由度:** rev1 は D₁(ℓ) と単一の家族レベル m を用いていたが、実際に定義できるのは
member ごとの D_{1,j}(ℓ) であり、27 member × 5 解像度から family verdict をどう決めるか（max_j、
member 平均、中央値、割合、per-member m_j の報告 …）が未指定だった。これは H₀ / H_dist の
最終判定を変えうる集約自由度である。

**A が登録した規則（§6.2.2、existential rule）と根拠:**

- member ごとに m_j = #{ℓ : D_{1,j}(ℓ) > τ_R}。member 判定 m_j = 0 / ≥2 / =1。
- **family H₀ ⟺ 全 member で m_j = 0；family H_dist ⟺ ∃j: m_j ≥ 2；さもなくば inconclusive。**
- 根拠は研究質問そのもの：情報階層は「E[V] が十分統計量か」を問い、**充足性は普遍的性質**である。
  したがって充足の family 主張は全 member で成り立つときのみ真であり、一つの member が
  （2 解像度以上で）平均を超える情報を要すれば反証される。多数決・平均・割合は「典型的 member」
  という別の質問（B が「科学的意味が異なる」と指摘したもの）に対応するため採らない。
- **τ_R 以外に新たなパラメータを導入しない**（割合も重みも member 選択もない）。
- **curve を見て選んでいない。**規則は充足性の意味論と、既計算の E[V], Var(V) のみから決めた。
  §6.2.6 のとおり、規則の定性的帰結（集中 member が H_dist を駆動しうること）はモーメントから
  透明に読めるが、距離 D_{k,j} 自体は凍結前に計算していない。

**K3/K4/K6 の機械的書き換え（§8）:** いずれも単一 m ではなく {m_j} → family verdict 集約規則に
依拠するよう改めた。K3 = family H₀、K4 = family inconclusive、K6 = embedding 間で family verdict
が相違。

### 11.4 本 RC(rev2) で B の最終確認を仰ぐ点

実質的な新規提案は §6.2.2 の集約規則 1 件のみ。確認を要するのは次の 2 点。

1. §6.2.2 の existential rule（family H₀ ⟺ 全 m_j=0；H_dist ⟺ ∃ m_j≥2）が、研究質問
   （E[V] の十分統計量性＝普遍性）に対応する member → family 規則として妥当か。
2. §6.2.6 の透明性開示（規則の帰結が既計算モーメントから定性的に読めること）が、
   「curve を見て規則を選んでいない」という要件と両立するか——すなわち、モーメントから
   予期される H_dist を開示することが outcome-based rule selection に当たらないという整理でよいか。

---

## 12. 変更履歴

- **v0.1–v0.3 DRAFT（2026-09-21）**: v1.0 §12 のとおり。凍結前。何も計算していない。
- **v1.0 FROZEN（2026-09-21）**: B 最終監査 PASS FOR FREEZE。何も計算していない。
- **v1.0 実装（2026-09-21）**: Haiku 4.5 版で T2–T5 がプレースホルダのまま PASS と記録され、
  P4/P6/P7/P8 が誤った対象を計算していた。`_INVALID_haiku45_2026-09-21/` に隔離。
  Opus 5 で再実装し T1–T5 を 27/27 通過（`PHASE_1_2_REPORT_v11.md`）。
- **P4 照会（2026-09-21）**: P4 が族内定数であることを発見。B 承認。
  E1b・E2 が登録形で計算不能。lesion 実行前。
- **解析評価（2026-09-21）**: redirect 介入自体の退化（20/27 member で η_½ = 1/2 厳密）、
  dead-end の十分統計量 (★) を導出。B PASS。lesion 実行前。
- **REPAIR SPEC v1.1（2026-09-21）**: B MINOR REVISION。3 点修正
  （closure の対称化、C3 の事後 matching 禁止、C6/E3 を control に留置）。
- **v1.1 DRAFT（2026-09-21）**: dead-end を primary に反転、情報階層 E[V] ⊂ {E[V],Var(V)} ⊂ P(V)、
  D 積分済み軌道推定量を primary、厳密 P(V) を oracle、removal-set を validation。
  凍結せず。何も outcome を計算していない。
- **v1.1 RELEASE CANDIDATE（2026-09-21）**: 本書。B REVIEW OF DRAFT
  "MINOR REVISION BEFORE FREEZE" の 4 点を反映（§11.2）。
  (1) E2 を三帰結に書き換え、D₂ を H_dist の必要条件から分類基準へ移した。
  (2) K4 を inconclusive / resolution-limited に書き換えた。
  (3) E1 を記述的診断に降格した。
  (4) T8 の 4 SE をグリッド点ごとの実装許容と明記した。
  他に実質的変更はない。凍結せず。何も outcome を計算していない。（この版を rev1 と呼ぶ。）
- **v1.1 RELEASE CANDIDATE rev2（2026-09-21）**: 本書。B FINAL AUDIT
  "ONE BLOCKING CLARIFICATION" に対応。§6.1 に member index j を明示し、§6.2 を
  **member → family の existential 集約規則**（family H₀ ⟺ 全 m_j=0；H_dist ⟺ ∃j m_j≥2；
  さもなくば inconclusive）に書き換えた。K3/K4/K6 を {m_j} 集約規則に合わせて機械的に更新。
  集約規則は充足性の意味論と既計算モーメントのみから決め、curve は計算していない（§11.3）。
  他に実質的変更はない。**凍結していない。何も outcome を計算していない。**
  registered removal seed 未消費。`tm_lesions_*.npz`、`tm_curves_*.npz` とも不存在。
- **v1.1 FROZEN（2026-09-21）**: B FINAL AUDIT **PASS FOR FREEZE**。RC rev2 に対し実質的変更なし
  （status header と本凍結宣言のみ追記）。凍結時点で outcome は何も計算していない。
  計算は §9 の順序に従い PI の明示的指示で開始する。

**すべての改訂は lesion outcome を見る前に行われた。改訂の根拠は経験的性能ではなく、
(★) の恒等式・§4.2 の導出・登録質問の意味論である。集約規則を含め、curve 値を見て
選んだ登録要素は一つもない。**
