# `.cursor/plans/` — 何がどのファイルか

このディレクトリは、**OrthoPaiR のリビジョン計画**と、その評価に用いる**ベンチマーク手順**の文書置き場です。実装コードではなく、方針・文献・ベンチマーク設計の文書です。

作業ブランチは **`revise_reorg`** を想定しています。

---

## 正本（まずこれを読む）

| 目的 | 読むファイル |
|------|----------------|
| **リビジョン計画の正本**（Phase 0–3 + 以降の publication） | [`revision_master_plan.md`](./revision_master_plan.md) |
| ベンチマーク手順の正典（Phase 0–8） | [`publication_grade_benchmark_protocol.md`](./publication_grade_benchmark_protocol.md) |
| ベンチマーク方法の文献調査（上記の背景） | [`orthology_benchmark_methods_2021_2026.md`](./orthology_benchmark_methods_2021_2026.md) |

**進捗・方針の「正」は常に [`revision_master_plan.md`](./revision_master_plan.md)。** 設計の詳細・経緯が必要なときだけ `archive/` を参照します。

```text
revision_master_plan.md ─┬─> Phase 0: ベンチマーク基盤の準備
（正本）                  ├─> Phase 1: 植物4ゲノム比較   = protocol Phase 1
                          ├─> Phase 2: アルゴリズム改修（Phase 1 結果駆動）
                          └─> Phase 3: 植物ベンチ評価 2.3 まで = protocol Phase 2 (§2.1–2.3)
                                   │
                                   └─(revision 完了後)─> publication 用 benchmarking
                                                          publication_grade_benchmark_protocol.md 全手順
```

---

## ファイル一覧（top-level）

| ファイル | 種類 | 一言 |
|----------|------|------|
| [`revision_master_plan.md`](./revision_master_plan.md) | 計画（正本） | リビジョン Phase 0–3 と、以降の publication 用 benchmarking への接続 |
| [`publication_grade_benchmark_protocol.md`](./publication_grade_benchmark_protocol.md) | 実施プロトコル | 査読論文に載せる評価の事前登録・階層・作業手順（Phase 0–8） |
| [`orthology_benchmark_methods_2021_2026.md`](./orthology_benchmark_methods_2021_2026.md) | 文献調査 | 2021–2026 のオーソログ評価方法論（QfO・OrthoBench・シミュレーション等） |
| [`README.md`](./README.md) | 索引 | 本ファイル |

---

## `archive/` — 正本に統合済みの設計・調査文書

以下は [`revision_master_plan.md`](./revision_master_plan.md) に内容を統合済みです。**設計の詳細・文献根拠・経緯**を確認する場合のみ参照してください。索引は [`archive/README.md`](./archive/README.md)。

| ファイル | 一言 |
|----------|------|
| [`archive/orthopair_algorithm_revision_plan.md`](./archive/orthopair_algorithm_revision_plan.md) | 03 OI/strictness・04 bridge/graph_qc の Phase 0–3 実装計画 |
| [`archive/orthopair_superiority_and_evaluation_plan.md`](./archive/orthopair_superiority_and_evaluation_plan.md) | 改修意図・現状・競合との差・現行ベスト数値 |
| [`archive/algorithm_from_preprints.md`](./archive/algorithm_from_preprints.md) | プレプリント手法 P1–P10 の適用検討・優先度 |
| [`archive/orthology_graph_heuristics_18901979.plan.md`](./archive/orthology_graph_heuristics_18901979.plan.md) | 狭義／広義・RBH グラフ近似を 03/04 のどこに載せるか |
| [`archive/reorg_giant_component_break_before_enumeration.md`](./archive/reorg_giant_component_break_before_enumeration.md) | 巨大連結成分を列挙前に壊す（初心者向け解説） |
| [`archive/revision_timeline_conflicts_and_direction.md`](./archive/revision_timeline_conflicts_and_direction.md) | 元計画 A と現行改修 B の時系列・衝突・統合方針 |

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

- **正本は [`revision_master_plan.md`](./revision_master_plan.md)。** 計画の変更はまず本書に反映する。
- **ベンチマーク手順の変更**は [`publication_grade_benchmark_protocol.md`](./publication_grade_benchmark_protocol.md) に反映する。
- **`archive/` の文書は原則追記しない。** 参照が必要になった内容は正本へ転載する。
- 実装の進捗とコードがずれたら、正本の該当 Phase を先に直す。
