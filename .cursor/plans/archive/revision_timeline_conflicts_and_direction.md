# 元計画と現行改修 — 時系列・衝突・今後の方針

**作成日:** 2026-09-04  
**ブランチ:** `aied`  
**目的:** もともとの Phase 計画と、実際に進んでいる reorg 回収サイクルを時系列で突き合わせ、衝突点を明示し、今後の方針案を固定する。

---

## 1. 二つの「改修計画」とは何か

| 呼称 | 主文書 | 北極星（一言） |
|------|--------|----------------|
| **元計画（A）** | [`orthopair_algorithm_revision_plan.md`](./orthopair_algorithm_revision_plan.md)（根拠: heuristics / preprints） | pseudoortholog 抑制・狭義/広義ラベル・多種グラフ QC（P4 / flooding）を **03→04→任意木** で段階実装 |
| **現行改修（B）** | [`REORG_REVISION_POLICY.md`](file:///home/ftom/01_wd/orthology/prior/cursor_dev/orthopair_revision_cycle/REORG_REVISION_POLICY.md) + cron cycle；整理は [`orthopair_superiority_and_evaluation_plan.md`](./orthopair_superiority_and_evaluation_plan.md) | link_all が示すグラフ上の BUSCO を、**選択的 reorg（A–G）**で回収。precision 不低下・ペア数爆発禁止 |

A は「何をアルゴリズムとして足すか」の設計書。  
B は「植物ベンチで OF/Sonic に追いつく」ための **運用駆動の実装ループ**。問題設定（Q1–Q3、二層構造）は共有するが、**優先順位と成功指標が分岐した**。

---

## 2. 時系列

```text
2026-05  A0  concepts
         orthology_graph_heuristics — 03/04 への載せ場所

2026-06-02  A1  literature
         algorithm_from_preprints — OI / P4 / flooding 優先

〜06初  A2  Phase plan 成立
         orthopair_algorithm_revision_plan — Phase 0–3

2026-06-03〜  B0  revision cycle 開始
         cycle 001: precision↑ recall↓ FAIL（過剰フィルタ）

2026-06-04  B1  北極星転換
         REORG_REVISION_POLICY 発動
         「OF/Sonic 級を選択的回収」が実装の正になる

2026-06-04〜08  B2  bridge 緩和・初回 PASS → PW ゲート達成

2026-06-08〜  B3  A–G 必須・飽和レバー禁止
         元計画 Phase 2 の「P4/flooding」より
         colourful / FN export / bridge path が主戦場

2026-06-16〜23  B4  プラトー + cycle 086+
         split-ID / hub orphan 等；BUSCO 横ばい帯

2026-07-04  統合メモ
         superiority_and_evaluation_plan
         （A の残タスク + B の現状を一文書に）

2026-07-18〜08-28  評価側の並行強化
         benchmark_methods → publication_grade_benchmark_protocol
         （論文用 D/V/H；開発用 cycle とは分離）

2026-09-03  計画文書の同期更新
         revision_plan / superiority / publication protocol

2026-09-04  現状認識
         実装コミットはほぼ B（A–G）
         A の YAML todos は Phase 2 in_progress のまま停滞
         巨大成分列挙コスト問題の解説メモも追加
```

### 2.1 元計画（A）で想定していた順序

| 順 | Phase | 内容 | 実装上の到達（2026-09 時点） |
|----|-------|------|------------------------------|
| 0 | ベンチ・baseline 固定 | 論文/ablation 用 | **未完了**（開発用 cycle eval が代替） |
| 1a | `shortest_path` のみ | 03 | **完了** |
| 1b | OI / strictness / flag_pseudo | 03 | **部分完了**（YAML は未更新） |
| 1c | IO / man / フィルタ API | 05 | **未完了** |
| 2 | bridge + graph_qc（弱辺・将来 P4） | 04 | bridge・弱辺 prune **あり**；P4 **未**。実装の大半は A–G |
| 2opt | 四分子 BLAST | 03/04 | **未着手** |
| 3 | 外部木エスカレーション | 別モジュール | **未着手** |

### 2.2 現行改修（B）で実際に起きた順序

| 順 | 出来事 | 成功指標 |
|----|--------|----------|
| 1 | pseudo フィルタ強化 → recall 崩壊 | FAIL |
| 2 | bridge 緩和・選択的回収 | 初回 PASS |
| 3 | PW 超えゲート固定 | commit 条件 |
| 4 | link_all 天井の発見 | 「落ちているだけ」仮説 |
| 5 | A–G 機構の積み上げ | gap 縮小＋precision 維持 |
| 6 | 飽和レバー禁止・プラトー | 微調整禁止；FN audit 先行 |
| 7 | superiority 文書で A 残件を再掲 | 優先はなお B（P0 = gap） |

---

## 3. 一致している点（衝突ではない）

次は両計画で同じなので、方針を戻す必要はない。

1. **二層構造** — ペア精度は 03、多種は 04。フル OrthoFinder 木は本体に入れない。
2. **`shortest_path` が唯一のペア経路** — 明示アンカー鎖は本線にしない。
3. **スコアは捨てずフラグ化** — `strictness` / `flag_pseudo` / OI。
4. **link_all 型の成分直積は禁止** — recall 天井の診断用のみ。
5. **Q1–Q3 の問題設定** — pseudoortholog、狭義/広義、多種一貫性。

---

## 4. コンフリクト（明示）

### C1. 成功指標の主軸

| | 元計画（A） | 現行（B） |
|--|------------|-----------|
| 主指標 | 1対1 precision↑、偽陽性（pseudo）↓、狭義ラベルの実用 | BUSCO recall を link_all/OF に近づけつつ precision≥PW、stretch≥Sonic |
| 緊張 | A は「偽 1対1 を落とす」方向が強い | B の初期過剰フィルタが FAIL → **落とすより回収**が主戦場に |

**衝突の本質:** A の Phase 1 精神（precision / 狭義）と、B の運用（recall 回収）は同じグラフ上でレバーが逆向きになりうる。

### C2. Phase 2 の中身

| | 元計画が書いた Phase 2 | 実際の Phase 2 相当作業 |
|--|------------------------|-------------------------|
| 中核 | P4 / near-cograph、弱辺 flooding、種系統スコア | colourful tie-break、FN 駆動 export、bridge path、残差第二パス（A–G） |
| graph_qc | 構造検定・ノイズ辺除去 | `prune_weak` 程度；巨大成分は割れない（別メモ） |

YAML の `phase2-graph-qc: in_progress` は **ラベル上は A の Phase 2**だが、コミット実体は **B の A–G**。文書と作業が同名で別物になっている。

### C3. 03 を触るタイミング

| | A | B（ポリシー） |
|--|---|---------------|
| 方針 | Phase 1 を最優先で 03 を先に固める | **reorg で埋まる限り `orthopair()` を変えない** |
| 結果 | 1b/1c が残ったまま | 04 に実装が集中 |

### C4. 評価の「正」

| | A / 論文プロトコル | B / cycle |
|--|-------------------|-----------|
| データ | D/V/Holdout、QfO、OrthoBench | 固定植物 4 ゲノム + BUSCO + PW ゲート |
| 危険 | cycle 数値を論文 primary に流用 | holdout 汚染・チューニング過適合 |

`publication_grade_benchmark_protocol.md` は「cycle gate だけでは論文主張に使わない」と明記。**衝突というより未統合**だが、数値の使い分けを誤ると方針衝突になる。

### C5. 計画文書の鮮度

- `orthopair_algorithm_revision_plan.md` の todos（1b pending 等）は **実装実態と不一致**。
- 進捗の正は現状 `orthopair_superiority_and_evaluation_plan.md` + git log（A–G）。
- 元計画を「今やっている計画」と読むと誤誘導（C2 と同根）。

### C6. 衝突ではないが未解決の技術ギャップ

- 巨大成分の colorful 列挙コスト（D1 煙テストで長時間固着）は、A の P4/flooding とも B の A–G とも別に、**列挙前に成分を壊す**必要が出ている（`reorg_giant_component_break_before_enumeration.md`）。これは両計画の前提（列挙が回る）を脅かす。

---

## 5. 今後の方針案

### 5.1 推奨：二層ロードマップ（統合案）

元計画を捨てず、**役割を分ける**。

```text
Track B（いま・短期）     Track A'（中期・製品/論文の差別化）
reorg 選択的回収 A–G      03/05 のラベル完成 + graph 構造 QC
+ 巨大成分の列挙前破壊     + Disputed / Tier C で「落としてよい偽陽性」を証明
        \                     /
         \                   /
          論文用評価 Track E
          publication_grade の D/V/H
          （cycle 数値は開発専用）
```

| Track | 期間感 | ゴール | 駆動文書 |
|-------|--------|--------|----------|
| **B** | 短期（継続） | link_all gap を precision 維持で縮小；ペア ≤2× PW；巨大成分を列挙可能にする | REORG_REVISION_POLICY、LINK_ALL_GAP、巨大成分メモ |
| **A'** | 中期 | Phase 1c 完了；`strict` を製品 API に；P4/`graph_qc` を **precision レバー**として再導入 | revision_plan（要更新）、preprints P1/P2 |
| **E** | 並行〜長期 | holdout・QfO・OrthoBench で査読可能な表 | publication_grade_benchmark_protocol |

### 5.2 衝突の解き方（ルール）

1. **Recall 回収（B）と precision 強化（A'）を同じコミットに混ぜない。**  
   - B 変更: commit gate = PW 超え + ペア数キャップ。  
   - A' の「落とす」変更: 別 PR/サイクル。評価は precision / Disputed / `strict` サブセットを必須。
2. **03 変更は B ポリシーどおり「reorg では埋められない FN」か「ラベルのみ（出力列・フィルタ）」に限定。**  
   - 1c（IO/man/strict フィルタ）はスコア分布を変えないので **いつでも可**。優先して潰してよい。
3. **元計画 YAML を実態に合わせて書き換える（文書コンフリクト C5 の解消）。**  
   - Phase 2 を「2a bridge+A–G（進行中）」「2b 構造 QC/P4（未着手）」「2c 列挙前成分破壊（新規）」に分割。
4. **論文の primary 数値は Track E のみ。** cycle / compare_tools 開発セットは Methods の development 節か supplementary。
5. **P4 / 弱辺強化は、B がプラトーかつ Disputed で「競合固有偽陽性」が示せたあと**に精度レバーとして入れる（A の精神を、B の天井が見えた後に回収）。

### 5.3 すぐやる順（提案）

| # | 作業 | Track | 理由 |
|---|------|-------|------|
| 1 | `revision_plan` の Phase 2 / todos を B 実態＋2b/2c に更新 | 文書 | 以後「今どの計画か」の誤読を止める |
| 2 | Phase 1c（`getOrthoPair`・man・strict フィルタ） | A' | スコア非破壊；Q2 のユーザー価値を先に出す |
| 3 | 巨大成分の列挙前破壊（prune 強化 / 分割） | B | 論文 D1 や大規模 reorg が回らないと評価不能 |
| 4 | FN audit 駆動の A–G 残り（ギャップが動くものだけ） | B | 飽和レバー禁止を維持 |
| 5 | Disputed pairs + `block_oi` 監査 | E/A' | 「Sonic 高 recall = 偽陽性混入」仮説の証拠 |
| 6 | P4 / graph_qc 拡張 | A' | precision stretch；B と別ゲート |
| 7 | publication holdout 1 回評価 | E | algorithm freeze 後 |

### 5.4 採らない案

| 案 | 却下理由 |
|----|----------|
| 元計画に完全回帰（A–G 停止、P4 最優先） | 植物ベンチで PW ゲートと gap 仮説が既に価値を出している；回帰は recall を捨てる |
| B のみ継続し A のラベル/QC を放棄 | 論文・製品の差別化（狭義 1対1、pseudo 説明）が消える |
| cycle 数値を論文 primary に昇格 | holdout 汚染；publication protocol と正面衝突 |
| link_all を deploy / 緩和した直積 | ポリシー禁止；ペア爆発 |

---

## 6. 一文での結論

**もともとの計画（A）は設計の骨格として生きているが、2026-06-04 以降の実装の正は reorg 選択的回収（B）に移った。**  
衝突の核は「偽陽性を落とす Phase 1 精神」と「グラフ上の真陽性を回収する A–G」のレバー方向である。今後は **B（回収＋巨大成分）→ A'（1c・strict・P4）→ E（論文 holdout）** に役割分担し、同じコミットで混ぜない。

---

## 7. 関連ドキュメント

| 文書 | 役割 |
|------|------|
| [orthopair_algorithm_revision_plan.md](./orthopair_algorithm_revision_plan.md) | 元計画 A（要: todos 更新） |
| [orthopair_superiority_and_evaluation_plan.md](./orthopair_superiority_and_evaluation_plan.md) | A+B 統合の現状ハブ |
| [REORG_REVISION_POLICY.md](file:///home/ftom/01_wd/orthology/prior/cursor_dev/orthopair_revision_cycle/REORG_REVISION_POLICY.md) | 現行 B の北極星 |
| [publication_grade_benchmark_protocol.md](./publication_grade_benchmark_protocol.md) | Track E |
| [reorg_giant_component_break_before_enumeration.md](./reorg_giant_component_break_before_enumeration.md) | 列挙前破壊（B の前提条件） |
| [README.md](./README.md) | 本ディレクトリ索引 |
