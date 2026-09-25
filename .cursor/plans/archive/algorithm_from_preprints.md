# プレプリント由来アルゴリズム案と OrthoPaiR への適用検討

**作成目的**: bioRxiv / arXiv 等で報告されているグラフ理論・オーソロジー解析の手法を整理し、OrthoPaiR で議論中の課題（狭義／広義オーソログ、pseudoortholog、RBH ネットワークによるツリー近似）に対して**何が使えそうか**、**精度向上が見込めるか**、**実装可能か**を検討する。

**関連プラン**: [orthology_graph_heuristics_18901979.plan.md](./orthology_graph_heuristics_18901979.plan.md)

**OrthoPaiR の現状（要約）**

| 段階 | ファイル | 要点 |
|------|----------|------|
| 01 | `R/01_functions_reformatFiles.R` | 入力正規化、最大 ~8999 ゲノム |
| 02 | `R/02_functions_rbh.R` | ペア／全ペア BLAST + RBH（`mutual_ci` 等） |
| 03 | `R/03_functions_orthology.R` | ペアごと（既定）: RBBH **アンカー** → C++ `rbh_shortest_path`（gene_index L1 距離）→ `.findShortestPath` → 1対1 化・分類（**オリジナルのシンテニー近似経路**）。オプションで明示アンカー鎖（`.findSyntenicOrtho`） |
| 04 | `R/04_functions_reorgOrthopairs.R` | 全ペア辺の `igraph`、**maximal colorful subgraph** 列挙、`orthopair_list.tsv` |
| 05 | `R/05_functions_IO.R` | 結果読出し |

---

## 1. 調査したプレプリント・近接文献（一覧）

| ID | 文献 | サーバー | グラフ／オーソロジーの要点 |
|----|------|----------|---------------------------|
| P1 | Lindeberg et al., *Orthology and Near-Cographs in the Context of Phylogenetic Networks* | [arXiv:2502.08746](https://arxiv.org/abs/2502.08746) | 真のオーソロジーグラフは cograph；ノイズ時は **near-cograph 編集**、level-1 ネットワークで説明可能なグラフの **線形時間認識** |
| P2 | Ma et al., *SOI: Robust identification of orthologous synteny with the Orthology Index* | [bioRxiv 2024.08.22.609065](https://www.biorxiv.org/content/10.1101/2024.08.22.609065v2) | シンテニーブロック内の **OI =（事前オーソログと判定された共線ペアの割合）**；OI≥0.6 で out-paralogy ブロック除去 |
| P3 | Klemm et al., *Proteinortho6: pseudo-reciprocal best alignment* | [Front. Bioinform. 2023](https://doi.org/10.3389/fbinf.2023.1322477)（査読済） | **疑似 RBH** で BLAST 半減；巨大 CC に **スペクトル分割 + 外れ値辺除去（flooding）** |
| P4 | Iwasaki & Cosentino, *SonicParanoid2* | [bioRxiv 2023.05.14.540736](https://www.biorxiv.org/content/10.1101/2023.05.14.540736v1) | AdaBoost で片方向 BLAST を省略；MCL；ドメイン Doc2Vec は別系統 |
| P5 | Emms et al., *OrthoFinder* (2025 update) | [bioRxiv 2025.07.15.664860](https://www.biorxiv.org/content/10.1101/2025.07.15.664860v1) | オーソグループ + **全遺伝子木 + duplication 解析**（参照・ベンチマーク用） |
| P6 | Altenhoff et al. 系 / FastOMA | [bioRxiv 2024.01.29.577392](https://www.biorxiv.org/content/10.1101/2024.01.29.577392v2) | 種系統誘導の HOG；**グラフ単体ではない**が多種一貫性の上限参照 |
| P7 | OrthoRefine | [BMC Bioinformatics 2024](https://doi.org/10.1186/s12859-024-05786-7) | 既存オーソグループを **シンテニー窓** で paralog 除去 |
| P8 | wQFM-DISCO | [bioRxiv 2023.12.05.570122](https://www.biorxiv.org/content/10.1101/2023.12.05.570122v1) | 多コピー遺伝子木を単コピー化して **四分子** 整合；種樹推定（オーソログペア直接ではない） |
| P9 | Dutilh / edge covering inparalog（古典） | [BMC Bioinformatics 2012 付録系](https://doi.org/10.1186/1471-2105-13-S19-S16) | 類似度グラフ上の **辺被覆** で in-paralog 近似（遺伝子木不要） |
| P10 | Synteny / phylogenomics (Brassicaceae) | [bioRxiv 2022.09.07.506897](https://www.biorxiv.org/content/10.1101/2022.09.07.506897v1) | シンテニー優先の **信頼できる 1対1**；遠縁でもブロック一貫性 |

※ medRxiv は本テーマ（ゲノムオーソロジー）ではほぼ該当なしのため省略。

---

## 2. OrthoPaiR の課題と文献アルゴリズムの対応

議論中の問題を次の3つに分解する。

1. **偽 1対1（pseudoortholog）**: 古い重複の別コピー同士が RBH で結ばれる。
2. **狭義 vs 広義**: 2種 MRCA レベルのオーソログ vs より古い重複由来の「包括的」相同。
3. **多種時のツリー代替**: 全遺伝子木なしで in/out-paralog に近い整合性。

以下、文献ごとに **アイデア → OrthoPaiR への載せ方 → 検証（精度・実装）** を記す。

---

### P1: Near-cograph / cograph 編集（arXiv 2502.08746）

**アルゴリズム要旨**

- 樹状進化の下ではオーソロジーグラフは **cograph**（誘導 P4 なし）。
- 推定グラフが cograph でないとき、(a) ノイズ、(b) ネットワーク進化、(c) **誤 RBH** のいずれか。
- **1 頂点削除で cograph になるグラフ = near-cograph**；level-1 ネットワークで説明可能なクラスを特徴づけ、**線形時間**で認識・ネットワーク構築。

**OrthoPaiR での使い方（アイデア）**

- **04 `reorgOrthopairs` 後**の連結成分（または疑わしい成分のみ）について:
  - 誘導 P4（a–b–c–d で a–c, b–d が辺）を列挙し、**4 遺伝子・4 種**のミスマッチ候補としてフラグ。
  - 辺削除コスト = `1 - mutual_ci` として **最小コストで cograph に近づける編集**（完全最適化は NP 困難な場合あり；ヒューリスティックで greedy 削除）。
- **03 ペア内**では P4 は情報不足になりやすい → 多種 04 向け。

**検証**

| 観点 | 評価 |
|------|------|
| 精度向上 | **中〜高（条件付き）**。誤 RBH が P4 を作る典型パターンと一致する場合、除去・再配分に有効。ただし論文の主眼は「グラフが進化史と整合するか」の**構造検定**であり、ペアワイズ精度のベンチマークではない。 |
| 実装可能性 | **中**。R では `igraph` で P4 検出・貪欲辺削除は容易。厳密な level-1 ネットワーク再構成までは過剰；**P4 フラグ + 弱辺削除**が現実的。 |
| 既存との関係 | 04 は既に **colorful subgraph**（種ごとに高々1遺伝子）を列挙しており、**「1種1遺伝子」の厳しいオーソグループ**方向と親和性が高い。near-cograph はその**前提グラフの妥当性チェック**層として追加できる。 |

---

### P2: Orthology Index（SOI, bioRxiv 2024）

**アルゴリズム要旨**

- シンテニーブロック内: **OI = n/m**（共線ペアのうち事前オーソログと呼ばれた割合）。
- OI≥0.6 で **orthologous synteny**；out-paralogy ブロックは低 OI。
- 下流で MCL により syntenic orthogroup（SOG）。

**OrthoPaiR での使い方（アイデア）**

- 03 で **ブロック単位 OI** を付与: 明示鎖経路（`_findSyntenyBlocks`）に加え、既定の `_findShortestPath`（アンカー + gene_index 距離）出力にも `strictness` / `flag_pseudo` を付ける。
- ペア TSV に `block_OI` 列；OI 低いブロック上の 1対1 は `class` を `broad` / `flag_outparalogy` に。
- 多種では 04 で同じ遺伝子ペアが複数ペア TSV に現れる → ブロック OI の **多種平均** を辺重みに。

**検証**

| 観点 | 評価 |
|------|------|
| 精度向上 | **高（植物・ポリプロイド）**。91 実データ + シミュレーションで orthology vs out-paralogy の分離が報告。遠縁2種でも **シンテニーが保たれていれば** 特に有効（ユーザー観察と一致）。 |
| 実装可能性 | **高**。OrthoPaiR は GFF・CDS 済み。ブロック定義と OrthoFinder 等の事前 OG は **自前の orthopair 結果で OI の分子・分母を置換可能**（SOI の思想のまま）。 |
| 注意 | OI は **事前オーソログの質に依存**（論文でも ~15% hidden paralog）。OrthoPaiR の mutual_ci フィルタと **併用**が前提。 |

---

### P3: Proteinortho6（疑似 RBH + スペクトルクラスタ + flooding）

**アルゴリズム要旨**

- **Pseudo-reciprocal**: 片方向 BLAST のみ実行し、逆方向 E-value を長さ・DB サイズから近似 → 計算半減。
- 巨大連結成分: **Fiedler ベクトルによる二分割**；無理なら **bitscore の外れ値辺を反復削除（flooding）** してから分割。
- 「全種に存在する core」仮定で適応的に分割継続。

**OrthoPaiR での使い方（アイデア）**

- **02**: 双方向 BLAST は既に最適化済みの可能性 → pseudo は **target_pair が疎なとき**のオプションに限定。
- **04**: `reorg` で巨大成分が `max_components` でスキップされる前に、**flooding 的辺削除**（`mutual_ci` が成分内中央値より Grubbs 的に低い辺）を試す。
- **03 では非推奨**: ペア内は CC が小さく、スペクトル分割の恩恵が小さい。

**検証**

| 観点 | 評価 |
|------|------|
| 精度向上 | **中**。巨大ハブクラスタの **チャンク化**には有効；個別 pseudoortholog 除去というより **誤マージ防止**。Proteinortho6 は QfO 等で BLAST 準拠と同等 precision を主張。 |
| 実装可能性 | **中〜高**。`igraph::cluster_leading_eigen` または軽量 flooding は R で可能。Lapack フル行列は大成分でメモリ爆発 → **辺リスト + サンプリング**必須。 |
| OrthoFinder 比較 | 計算削減は **BLAST 段階**；OrthoPaiR のボトルネックが BLAST なら有益、03–04 だけなら効果小。 |

---

### P4: SonicParanoid2（ML で BLAST 省略 + MCL）

**アルゴリズム要旨**

- AdaBoost で **速い方向の inter-proteome alignment を先に実行**し、essential subset だけ逆方向。
- MCL で OG；ドメイン Doc2Vec は融合遺伝子向け。

**OrthoPaiR での使い方（アイデア）**

- 02 の **ジョブ順序最適化**（大きいゲノムを query にする等）の理論的根拠として参考。
- **本質的な狭義／広義分離には直接寄与しにくい**（OG 全体の recall 向上が主）。

**検証**

| 観点 | 評価 |
|------|------|
| 精度向上 | **低〜中**（ユーザーの pseudoortholog 問題に対して）。QfO で Pareto 最適だが、**1対1 遠縁ペア**の評価ではない。 |
| 実装可能性 | **低**（AdaBoost モデル再学習・パイプライン分岐）。既存 RBH パイプラインとの統合コスト大。 |
| 推奨 | **スケール時の BLAST 削減**のみ検討；オーソログ「質」の主戦場は 03/04 + P2/P1。 |

---

### P5: OrthoFinder 2025（遺伝子木 + duplication）

**アルゴリズム要旨**

- オーソグループ → 遺伝子木 → duplication マッピング；精度 +7% 等。

**OrthoPaiR での使い方**

- **ゴールド標準ベンチマーク**および **段階的エスカレーションの行き先**（疑わしいペアだけ OrthoFinder サブセット）。
- OrthoPaiR 本体に木推定を載せると **設計思想が OrthoFinder 化** → 別モジュール推奨。

**検証**

| 観点 | 評価 |
|------|------|
| 精度向上 | **高（上限）** |
| 実装可能性 | **別パイプライン**として現実的；組み込みは **低**（コスト・依存関係） |

---

### P6: FastOMA（種系統 + HOG）

**参考**: 千ゲノム規模では **種系統上の葉→根トラバーサル**で HOG を解く。グラフヒューリスティックの代替ではなく、**taxonomy-guided の経路一貫性**の実装例。

**OrthoPaiR**: ユーザーが種系統を渡せば、04 で **辺 (g1,g2) が種系統の隣接／祖先関係と矛盾しないか**をスコア化可能。**実装は中**、精度は系統の正しさに依存。

---

### P7: OrthoRefine（シンテニー窓での paralog 除去）

**アルゴリズム要旨**

- 既存 OG の各遺伝子について近傍ウィンドウ内の **シンテニー保存**を見て paralog を落とす。

**OrthoPaiR での使い方**

- 03 `_filterOrthopair` の拡張、または 04 で **同一 OG 候補内の複数コピー**をウィンドウで間引く。
- P2（OI）と **ほぼ同系統**；OI はブロック単位、OrthoRefine は遺伝子単位。

**検証**

| 観点 | 評価 |
|------|------|
| 精度向上 | **中〜高**（細菌・真核データで specificity 向上） |
| 実装可能性 | **高**（GFF 上の look-around） |

---

### P8: 四分子・DISCO（系統ゲノミクス系）

**アルゴリズム要旨**

- 遺伝子木の多コピーを分解し **四分子分布**を種樹と整合。

**OrthoPaiR での使い方**

- BLAST 上位ヒットから **4 種 (A,B,C,O)** の「AB|CD vs AC|BD」的整合を **スコア化**（木を立てない）。
- 03 で `mutual_ci` が低いペアに限定すれば **計算は O(疑わしいペア数)**。

**検証**

| 観点 | 評価 |
|------|------|
| 精度向上 | **中**（外群 1 種で古い重複のすれ違い検出に有効） |
| 実装可能性 | **中**（第三・第四ゲノムの RBH が必要） |

---

### P9: 辺被覆による in-paralog（グラフ理論）

**アルゴリズム要旨**

- 同一種内の in-paralog を **類似度グラフの被覆**で選び、種間オーソログ代表とマッチング。

**OrthoPaiR での使い方**

- 03 の `_pickBestPair` 前に、**種内 RBH クラスタ**から代表転写産物を1本に（02 で種内 BLAST があれば）。
- 多種 04 では **colorful subgraph** と併用し、同一種内の複数頂点を被覆で1つに。

**検証**

| 観点 | 評価 |
|------|------|
| 精度向上 | **中**（1対多・多対多の整理；偽1対1の一部に効く） |
| 実装可能性 | **中**（種内 hit が無いと効果限定；02 拡張が必要な場合あり） |

---

### P10: シンテニー優先 phylogenomics

**アルゴリズム要旨**

- 遠縁でも **共線ブロック上の 1対1** を系統解析に採用。

**OrthoPaiR での使い方**

- **既定の `_findShortestPath` + `rbh_shortest_path` が、すでにシンテニー（アンカー RBH + 共線 chr 上の gene_index 順序）を利用するオリジナル経路**である点を区別する。
- 追加オプションとして **`_findSyntenicOrtho`（明示アンカー鎖）** の比較・`hybrid`（高信頼のみ明示鎖、残りオリジナル）が「既存資産の活用」。

**検証**

| 観点 | 評価 |
|------|------|
| 精度向上 | **高（ユーザー実績ベース）** — 主にオリジナル経路の強化（OI ラベル）と、必要時のみ明示鎖 |
| 実装可能性 | **高**（オリジナルは稼働中；明示鎖はオプション） |

---

## 3. 統合評価マトリクス

**問題**: (A) pseudoortholog 除去 (B) 狭義／広義ラベル (C) 多種・ツリー代替

| 案 | 主文献 | A | B | C | 実装先 | 実装難度 | 精度見込み |
|----|--------|---|---|---|--------|----------|------------|
| OI / ブロックシンテニー | P2, P7, P10 | ◎ | ◎ | ○ | 03 (+04) | 低〜中 | 高 |
| オリジナル経路の OI/ラベル + 明示鎖オプション | P10, archive/05 | ◎ | ◎ | ○ | 03 | 低 | 高（2種） |
| P4 / near-cograph 辺編集 | P1 | ○ | ○ | ◎ | 04 | 中 | 中〜高 |
| Flooding / スペクトル分割 | P3 | △ | △ | ◎ | 04 | 中 | 中（マージ誤り） |
| 四分子 BLAST 制約 | P8 | ○ | ◎ | ◎ | 03 | 中 | 中 |
| 種系統パス一貫性 | P6 | ○ | ◎ | ◎ | 04 | 中 | 中（系統依存） |
| 辺被覆 in-paralog | P9 | ○ | ○ | ○ | 02–03 | 中 | 中 |
| ML BLAST 省略 | P4 | △ | △ | △ | 02 | 高 | 低（本問題） |
| 全遺伝子木・調停 | P5 | ◎ | ◎ | ◎ | 外部 | 高 | 高（上限） |

◎=強く有効 ○=補助 △=副次的

---

## 4. 推奨ロードマップ（実装優先度）

### Phase 1 — 低コスト・高精度が見込める（03 中心）

1. **既定 `_findShortestPath`（オリジナル・シンテニー近似）を維持**しつつ、ペア内 **Orthology Index / `strictness`**（P2）を列追加。
2. **オプション**: `_findSyntenicOrtho` / `_findSyntenyBlocks`（明示鎖）と `hybrid`（P10, archive/05 比較用）。
3. **`class` の二段化**: `1to1_strict` / `1to1_broad` / `flag_pseudo`（閾値はベンチマークで校正）。

### Phase 2 — 多種グラフ（04 中心）

4. **疑わしい成分のみ P4 検出 + 弱 `mutual_ci` 辺の greedy 削除**（P1）。
5. **巨大成分への flooding**（P3）を `max_components` スキップの前処理に。
6. **種系統（任意入力）による辺妥当性スコア**（P6, P8 の軽量版）。

### Phase 3 — 検証とエスカレーション

7. 小さな **QfO / シミュレーション / OrthoFinder HOG** で Phase 1–2 の precision–recall。
8. `flag_pseudo == TRUE` のみ **FastTree + 簡易調停** または OrthoFinder サブセット（P5）。

---

## 5. 検証方法（実装前に設計すべきベンチマーク）

プレプリントはタスクが異なるため、OrthoPaiR 向けには次を推奨する。

1. **ペアワイズ**: 遠縁2種で手動／文献カーブーションした「真の 1対1」セット（数十〜数百遺伝子）。
2. **多種**: OrthoFinder または SOI パイプライン出力を **参照 OG** とし、OrthoPaiR の `orthopair_list.tsv` 行の precision/recall。
3. **Ablation**: (i) RBH のみ (ii) +オリジナル shortest_path (iii) +OI/strictness 列 (iv) +明示鎖 or hybrid (v) +04 graph edit。
4. **指標**: 1対1 precision（偽陽性＝pseudoortholog）、recall、および `mutual_ci` との相関。

**期待される結果（仮説）**

- Phase 1 だけで **偽 1対1 precision が有意に向上**（ユーザーの既報と P2/P10 が支持）。
- P1/P3 は **巨大クラスタの誤結合**に効くが、2種のみ解析では効果が小さい場合がある。
- 木ベース（P5）を超える recall は **期待しない**；狭義 precision の向上が現実的目標。

---

## 6. 結論

- プレプリント群のうち、OrthoPaiR の**既存アーキテクチャに最も自然に載る**のは **P2（OI）・P10（シンテニー）・P1（グラフ構造検定）・P3（巨大成分の辺整理）**。
- **ツリー調停に近い狭義／広義の区別**は、フル木なしでは **OI + 四分子／種系統制約 + near-cograph 編集**の組み合わせで **近似**するのが文献とも整合的。
- **実装可能性が高く精度向上が見込める**のは Phase 1（03 **オリジナル経路のラベル強化**＋必要時の明示鎖）。**04 へのグラフ理論**は中程度の実装コストで多種・偽マージに効く。
- **SonicParanoid2 型 ML** や **OrthoFinder 級の全木**は、本問題の主解としては優先度低（スケール／別製品線）。

---

## 7. 参考文献（URL）

- P1: https://arxiv.org/abs/2502.08746  
- P2: https://www.biorxiv.org/content/10.1101/2024.08.22.609065v2  
- P3: https://doi.org/10.3389/fbinf.2023.1322477  
- P4: https://www.biorxiv.org/content/10.1101/2023.05.14.540736v1  
- P5: https://www.biorxiv.org/content/10.1101/2025.07.15.664860v1  
- P6: https://www.biorxiv.org/content/10.1101/2024.01.29.577392v2  
- P7: https://doi.org/10.1186/s12859-024-05786-7  
- P8: https://www.biorxiv.org/content/10.1101/2023.12.05.570122v1  
- P10: https://www.biorxiv.org/content/10.1101/2022.09.07.506897v1  
