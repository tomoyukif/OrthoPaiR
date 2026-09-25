# `.cursor/plans/` — 何がどのファイルか

このディレクトリは、**OrthoPaiR アルゴリズム改修**（狭義／広義オーソログ、pseudoortholog 抑制、03/04 改修、論文用評価）の計画・調査メモ置き場です。実装コードではなく、方針・文献・ベンチマーク設計の文書です。

作業ブランチは **`revise_reorg`** を想定しています。ファイル名は他文書・ワークスペースルールから参照されているため、整理は README による索引とし、リネームはしません。

---

## まずどれを読むか

| 目的 | 読むファイル |
|------|----------------|
| 元計画 vs 現行改修の時系列・衝突・今後方針 | [`revision_timeline_conflicts_and_direction.md`](./revision_timeline_conflicts_and_direction.md) |
| 問題設定とパイプラインの置き場所 | [`orthology_graph_heuristics_18901979.plan.md`](./orthology_graph_heuristics_18901979.plan.md) |
| 文献から何を実装に載せるか | [`algorithm_from_preprints.md`](./algorithm_from_preprints.md) |
| Phase 0–3 の実装計画 | [`orthopair_algorithm_revision_plan.md`](./orthopair_algorithm_revision_plan.md) |
| 改修の意図・現状・競合との差 | [`orthopair_superiority_and_evaluation_plan.md`](./orthopair_superiority_and_evaluation_plan.md) |
| 論文に足る評価の手順 | [`publication_grade_benchmark_protocol.md`](./publication_grade_benchmark_protocol.md) |
| ベンチマーク方法の文献調査 | [`orthology_benchmark_methods_2021_2026.md`](./orthology_benchmark_methods_2021_2026.md) |
| 巨大成分で reorg が終わらない理由と壊し方（初心者向け） | [`reorg_giant_component_break_before_enumeration.md`](./reorg_giant_component_break_before_enumeration.md) |

推奨順（概念 → 実装 → 評価）:

1. `orthology_graph_heuristics_18901979.plan.md`
2. `algorithm_from_preprints.md`
3. `orthopair_algorithm_revision_plan.md`
4. `orthopair_superiority_and_evaluation_plan.md`
5. `orthology_benchmark_methods_2021_2026.md`
6. `publication_grade_benchmark_protocol.md`

```text
概念・置き場所          文献適用               実装計画
heuristics.plan.md  →  algorithm_from_preprints.md  →  algorithm_revision_plan.md
                                                         ↓
文献調査（評価方法）                                    優位性・現状整理
benchmark_methods_2021_2026.md  →  superiority_and_evaluation_plan.md
                                                         ↓
                                                   論文用ベンチ手順
                                                   publication_grade_benchmark_protocol.md
```

---

## ファイル一覧

| ファイル | 種類 | 一言 |
|----------|------|------|
| [`revision_timeline_conflicts_and_direction.md`](./revision_timeline_conflicts_and_direction.md) | 方針整理 | 元計画 A と現行 B の時系列・衝突・今後の Track 分担 |
| [`orthology_graph_heuristics_18901979.plan.md`](./orthology_graph_heuristics_18901979.plan.md) | 概念プラン（Cursor YAML） | 01–05 パイプライン上で、狭義／広義と RBH グラフ近似をどこに載せるか |
| [`algorithm_from_preprints.md`](./algorithm_from_preprints.md) | 文献調査 | プレプリント手法を OrthoPaiR 課題に対応づけ、実装優先度を付ける |
| [`orthopair_algorithm_revision_plan.md`](./orthopair_algorithm_revision_plan.md) | 実装計画（Cursor YAML） | Phase 0–3：03 OI/strictness、04 bridge/graph_qc |
| [`orthopair_superiority_and_evaluation_plan.md`](./orthopair_superiority_and_evaluation_plan.md) | 現状整理＋評価方針 | 改修の意図、実装状況、競合との差、残作業 |
| [`orthology_benchmark_methods_2021_2026.md`](./orthology_benchmark_methods_2021_2026.md) | 文献調査 | 2021–2026 のオーソログ評価方法論（QfO、OrthoBench、シミュレーション等） |
| [`publication_grade_benchmark_protocol.md`](./publication_grade_benchmark_protocol.md) | 実施プロトコル | 査読論文に載せる評価の事前登録・階層・作業手順 |
| [`reorg_giant_component_break_before_enumeration.md`](./reorg_giant_component_break_before_enumeration.md) | 解説メモ（初心者向け） | 巨大連結成分で colourful 列挙が終わらない理由と、列挙前に壊す手順 |

`.plan.md` は Cursor のプラン形式（先頭 YAML に `todos`）です。通常の `.md` は調査・方針メモです。

---

## 各ファイルの内容

### `revision_timeline_conflicts_and_direction.md`

**作成:** 2026-09-04  
**役割:** 元計画（Phase 0–3）と現行 reorg 回収サイクルの時系列比較、コンフリクト明示、今後の Track B / A' / E 分担案。

### `orthology_graph_heuristics_18901979.plan.md`

**作成:** 2026-05-21 頃  
**役割:** 改修の出発点。OrthoPaiR の 01–05 パイプラインを前提に、狭義／広義オーソログと「RBH ネットワークでツリーを近似する」話を、**03（ペア内）と 04（多種グラフ）のどこに載せるか**を整理する。

含むもの:

- コア5ファイル（reformat → RBH → orthopair → reorg → IO）の対応表
- 既定経路 `shortest_path`（RBBH アンカー + gene_index L1）がシンテニー近似であること
- 03 の行メタデータ vs 04 の後段スコア vs 疑わしいものだけ木へ送る二段処理

後続の実装計画はこれを前提にしている。

### `algorithm_from_preprints.md`

**作成:** 2026-06-02 頃  
**役割:** bioRxiv / arXiv 等の手法を、OrthoPaiR の3課題（pseudoortholog、狭義 vs 広義、多種時のツリー代替）に対応づける。

含むもの:

- 文献一覧 P1–P10（near-cograph、SOI/OI、Proteinortho6 flooding、SonicParanoid2、OrthoFinder、FastOMA 等）
- 各手法の「OrthoPaiR への載せ方」と実装可能性
- 推奨ロードマップ（Phase 1: OI/シンテニー、Phase 2: 04 グラフ QC、Phase 3: 任意の木調停）

`orthopair_algorithm_revision_plan.md` の文献根拠。

### `orthopair_algorithm_revision_plan.md`

**最終更新:** 2026-09-03 頃  
**役割:** 上記2文書を実装タスクに落とした **Phase 0–3 計画**。ブランチ `revise_reorg`。

| Phase | 内容 | メモ上の状態（YAML todos） |
|-------|------|---------------------------|
| 0 | ベンチマーク設計と現行 baseline | pending |
| 1a | `orthopair()` を `shortest_path` のみに | completed |
| 1b | ブロック OI・strictness・`flag_pseudo` | pending（後続文書では部分完了） |
| 1c | `getOrthoPair`・man・後方互換 | pending |
| 2 | 04 の bridge・弱辺 flooding（graph_qc） | in_progress |
| 2 任意 | 四分子 BLAST | pending |
| 3 任意 | `flag_pseudo` の外部木フック | pending |

スコープ外: OrthoFinder 級の全遺伝子木、Sonic 型 ML による BLAST 省略、DIAMOND 全面置換。

YAML の todo 状態と、`orthopair_superiority_and_evaluation_plan.md` の「実装状況」表は更新タイミングがずれうる。進捗の正は後者を優先して確認する。

### `orthopair_superiority_and_evaluation_plan.md`

**作成:** 2026-07-04（2026-09-03 頃更新）  
**役割:** 「なぜ改修するか」「今どこまでできているか」「OrthoFinder / Sonic に対する優位をどう示すか」を一箇所にまとめたハブ文書。

含むもの:

- 設計原則（二層構造、スコアは捨てずフラグ化、link_all 型直積禁止、選択的 reorg）
- reorg 北極星（recall は OF 天井に近づける、precision は落とさない、ペア数 ≤ 2× PW）
- Phase 1/2 の実装状況と reorg 機構 A–G
- 植物ベンチの現行ベストと競合比較の位置づけ
- 残りのアルゴリズム改修とツール比較の方法論

実装計画と論文用プロトコルの橋渡し。reorg の北極星の詳細はリポジトリ外の  
`orthology/prior/cursor_dev/orthopair_revision_cycle/REORG_REVISION_POLICY.md`。

### `orthology_benchmark_methods_2021_2026.md`

**作成:** 2026-07-18（2026-07-28 頃更新）  
**役割:** 「精度」を何で測るかの**方法論調査**。実装計画ではなく、評価設計の文献根拠。

含むもの:

- QfO、VGNC、SwissTree、Open OrthoBench、STD、GO/EC、FAS、シミュレーション、シンテニー/WGD、ツール間不一致、計算性能
- OrthoPaiR 改修前に推奨する最小構成・植物特化セット・ablation
- BUSCO ペアの precision/recall だけでは論文主張に足りない、という結論

`publication_grade_benchmark_protocol.md` の背景資料。

### `publication_grade_benchmark_protocol.md`

**作成:** 2026-08-28（2026-09-03 頃更新）  
**役割:** 査読論文で主張できる水準まで評価するための **方針と実施手順**。

含むもの:

- Development / Validation / Holdout の分割と、holdout は freeze 後に1回だけ
- Primary / secondary endpoints、比較ツールと版固定
- 論文 Main vs Supplementary の階層（植物 holdout、OrthoBench、runtime、QfO 等）
- 成果物ディレクトリ案、フェーズ別作業、投稿前品質ゲート

開発用 `orthopair_revision_cycle` の commit gate だけでは論文主張に使わない、と明記している。

### `reorg_giant_component_break_before_enumeration.md`

**作成:** 2026-09-04  
**役割:** D1 煙テストで `reorg` が CPU 張り付きに見えた事象をきっかけに、**巨大連結成分・colorful 列挙・列挙前に壊す**を初心者向けに説明する。

含むもの:

- 連結成分 / colorful subgraph / Bron–Kerbosch の直感的定義
- `prune_weak` だけでは巨大成分が割れない理由と、D1 の数値例
- 本命（1to1 backbone で列挙成分を切り直す）、保険（反復 flooding）、例外（`max_components`）
- link_all 直積などやってはいけないこと
- 実装の優先順と API 案

Phase 2 の graph_qc / flooding 議論の入門。実装タスクの正本は `orthopair_algorithm_revision_plan.md`。

---

## このディレクトリの外（よく参照される関連文書）

| 場所 | 内容 |
|------|------|
| `orthology/prior/cursor_dev/orthopair_revision_cycle/REORG_REVISION_POLICY.md` | reorg 改修の北極星・コミットゲート |
| `orthology/prior/cursor_dev/compare_tools/compare_tools_plan.md` | 多層証拠によるツール差分評価 |
| `orthology/prior/cursor_dev/validate_ortho/evaluation_method.md` | Tier C（遺伝子樹）検証 |

パスは OrthoPaiR リポジトリではなく、隣の `orthology` ワークスペース側です。

---

## メンテ方針

- **追加するとき:** この README の一覧表と「各ファイルの内容」に1段落足す。
- **役割が重なったら:** リネームせず、ハブ側（現状は `orthopair_superiority_and_evaluation_plan.md`）から「正」をどちらにするか書いて相互リンクする。
- **実装の進捗:** コードと YAML todos がずれたら、`orthopair_superiority_and_evaluation_plan.md` の実装状況表を先に直す。
