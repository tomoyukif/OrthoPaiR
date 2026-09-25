# `.cursor/plans/archive/` — アーカイブ済み計画・調査文書

ここにある文書は、**2026-09-25 のリビジョン計画の再編**で正本 [`../revision_master_plan.md`](../revision_master_plan.md) に内容を統合したため、top-level から移動したものです。

**これらは「今の計画」ではありません。** 進捗・方針の正は常に [`../revision_master_plan.md`](../revision_master_plan.md) を参照してください。ここは **設計の詳細・文献根拠・経緯**を辿るための参照用です。

---

## なぜアーカイブしたか / 正本のどこに吸収されたか

| ファイル | 種類 | 正本での吸収先 |
|----------|------|----------------|
| [`orthopair_algorithm_revision_plan.md`](./orthopair_algorithm_revision_plan.md) | 実装計画（Phase 0–3） | 正本 §4（Phase 2 アルゴリズム改修の詳細タスク・API 案・列定義） |
| [`orthopair_superiority_and_evaluation_plan.md`](./orthopair_superiority_and_evaluation_plan.md) | 現状整理＋評価方針ハブ | 正本 §4.2（北極星）・§3.3（現行ベスト数値）・競合位置づけ |
| [`algorithm_from_preprints.md`](./algorithm_from_preprints.md) | 文献調査（P1–P10） | 正本 §4.4（プレプリント根拠の要約表） |
| [`orthology_graph_heuristics_18901979.plan.md`](./orthology_graph_heuristics_18901979.plan.md) | 概念プラン | 正本 §4.1–§4.2（二層構造・問題設定） |
| [`reorg_giant_component_break_before_enumeration.md`](./reorg_giant_component_break_before_enumeration.md) | 解説メモ | 正本 §4.5（巨大成分の列挙前破壊） |
| [`revision_timeline_conflicts_and_direction.md`](./revision_timeline_conflicts_and_direction.md) | 方針整理（A/B 衝突） | 正本 §4.6（改修の進め方・衝突回避ルール） |

---

## 注意

- これらの文書内の相対リンク（`./` や `../../`）は、移動前（`.cursor/plans/` 直下）を前提にしているため、一部リンクが切れている場合があります。
- 記載の作業ブランチは `aied` 想定の記述が残りますが、現行の作業ブランチは `revise_reorg` です（正本を参照）。
