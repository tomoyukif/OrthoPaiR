# オーソログ解析ツールのベンチマーク方法：2021–2026 年の文献調査

**作成日:** 2026-07-18  
**対象期間:** 原則 2021–2026 年（方法論の基礎として 2016・2020 年の標準文献も補足）  
**目的:** OrthoPaiR の追加改修に先立ち、何を「精度」と定義し、どのベンチマークを採用すべきか整理する。

---

## 1. 結論

オーソログ推定には、全遺伝子・全系統を覆う完全な ground truth がない。そのため、近年の標準的な評価は、次の独立性の異なる評価を組み合わせる。

1. **共通入力上の標準ベンチマーク** — Quest for Orthologs（QfO）
2. **専門家キュレーション済み正解集合** — VGNC、SwissTree、TreeFam-A、Open OrthoBench
3. **既知の種系統との整合性** — Species Tree Discordance
4. **機能・タンパク質構造特徴の保存** — GO、EC、Feature Architecture Similarity（FAS）
5. **シミュレーションによる既知の真値** — ALF、Zombi、DTL/ILS reconciliation sims 等
6. **シンテニー・倍数性・遺伝子順序の整合性** — 植物・WGD 系統に特に重要
7. **ツール間不一致と頑健性** — consensus、Signal Jaccard Index、入力摂動
8. **下流解析への影響** — 種系統、遺伝子 gain/loss、phylostratigraphy
9. **計算性能** — wall time、CPU、peak RSS、I/O、スケーリング、完走率

したがって、OrthoPaiR の優位性を示すベンチマークは、BUSCO ペアの precision / recall だけでは不十分である。最低でも **QfO + Open OrthoBench + 植物 WGD/シンテニー + シミュレーション + 計算性能**を採用し、開発用データと最終評価用 holdout を分ける必要がある。STD の RF は precision ではなく **topological concordance proxy** として扱い、event-level の真値が必要な場合は curated trees か simulation/reconciliation を併用する。

---

## 2. 近年の主要文献

### 2.1 調査対象文献

| 年 | 文献 | ベンチマーク上の主な貢献 |
|---:|------|--------------------------|
| 2021 | Natsidis et al., *Systematic errors in orthology inference and their effects on evolutionary analyses*, iScience 24:102110. [DOI](https://doi.org/10.1016/j.isci.2021.102110) | ALF シミュレーションで真の ortholog を既知にし、進化速度依存の系統的誤りと下流解析への影響を評価 |
| 2021 | Conover et al., *pSONIC: Ploidy-aware Syntenic Orthologous Networks Identified via Collinearity*, G3 11:jkab170. [DOI](https://doi.org/10.1093/g3journal/jkab170) | 植物・倍数体で、single-copy gene tree の種系統一致率と、倍数体をサブゲノム分割した場合の出力一致率を評価 |
| 2022 | Nevers et al., *The Quest for Orthologs orthology benchmark service in 2022*, Nucleic Acids Research 50:W623–W632. [DOI](https://doi.org/10.1093/nar/gkac330) | QfO 共通プロテオーム、OpenEBench、VGNC の人手 curated ortholog benchmark を追加 |
| 2023 | Dosch et al., *FAS: assessing the similarity between proteins using multi-layered feature architectures*, Bioinformatics 39:btad226. [DOI](https://doi.org/10.1093/bioinformatics/btad226) | ドメイン・膜貫通領域・低複雑性領域などの feature architecture 保存度を 0–1 で評価 |
| 2024 | Altenhoff et al., *New developments for the Quest for Orthologs benchmark service*, NAR Genomics and Bioinformatics 6:lqae167. [DOI](https://doi.org/10.1093/nargab/lqae167) | FAS を QfO の新しい標準ベンチマークとして追加し、18 手法を共通入力上で比較 |
| 2024 | Cosentino et al., *SonicParanoid2: fast, accurate, and comprehensive orthology inference with machine learning and language models*, Genome Biology 25:195. [DOI](https://doi.org/10.1186/s13059-024-03298-4) | QfO の複数精度試験と、QfO・200 eukaryotes・2,000 MAGs による速度/メモリ/スケーラビリティ評価 |
| 2024 | Ludwig and Mrázek, *OrthoRefine: automated enhancement of prior ortholog identification via synteny*, BMC Bioinformatics 25:163. [DOI](https://doi.org/10.1186/s12859-024-05786-7) | 既存 orthogroup をシンテニーで精緻化し、近傍保存を独立証拠として利用 |
| 2025 | Majidian et al., *Orthology inference at scale with FastOMA*, Nature Methods. [DOI](https://doi.org/10.1038/s41592-024-02552-8) | QfO 11 試験、Bilateria orthogroup benchmark、2,086 eukaryotes の大規模性能を併用 |
| 2025 | Yang et al., *A metric and its derived protein network for evaluation of ortholog database inconsistency*, BMC Bioinformatics 26:6. [DOI](https://doi.org/10.1186/s12859-024-06023-x) | genome-context に基づく Signal Jaccard Index（SJI）とネットワーク degree centrality による不一致・信頼性評価 |
| 2026 | Emms et al., *OrthoFinder: improved phylogenetic orthology inference with enhanced accuracy and scalability*, Nature Methods. [DOI](https://doi.org/10.1038/s41592-026-03126-6) | QfO と expert-curated OrthoBench、orthogroup の fission/fusion、計算量を組み合わせた最新比較 |
| 2022 | Morel et al., *SpeciesRax…*, Mol Biol Evol. [DOI](https://doi.org/10.1093/molbev/msab365) | DTL シミュレーションで種系統再構成精度を評価（reconciliation 系の event-level 真値） |
| 2024 | Morel et al., *AleRax…*, Bioinformatics. [DOI](https://doi.org/10.1093/bioinformatics/btae162) | gene-tree 不確実性を含む DTL reconciliation 精度評価 |
| 2023– | Mishra et al., *reconcILS* / dupcoal, bioRxiv. [DOI](https://doi.org/10.1101/2023.11.03.565544) | ILS + duplication/loss 同時シミュレーションで、STD が偽陽性とみなす discordance を event 単位で分離 |

### 2.2 期間外だが外せない基礎文献

| 年 | 文献 | 理由 |
|---:|------|------|
| 2016 | Altenhoff et al., *Standardized benchmarking in the quest for orthologs*, Nature Methods 13:425–430. [DOI](https://doi.org/10.1038/nmeth.3830) | QfO の species-tree、reference-tree、functional benchmark の定義 |
| 2020 | Emms and Kelly, *Benchmarking Orthogroup Inference Accuracy: Revisiting Orthobench*, Genome Biology and Evolution 12:2258–2266. [DOI](https://doi.org/10.1093/gbe/evaa211) | 現在も利用される Open OrthoBench。70 個の expert-curated Bilaterian RefOG と評価スクリプトを公開 |

Open OrthoBench は 2020 年公開で「過去 5 年」からわずかに外れるが、FastOMA（2025）と OrthoFinder（2026）が現在も採用しているため、本調査に含める。

---

## 3. ベンチマーク方法の分類

## 3.1 Quest for Orthologs（QfO）標準ベンチマーク

### 方法

QfO Reference Proteomes の**同一配列・同一 ID**を全ツールに入力し、予測を OrthoXML または protein pair TSV として OpenEBench に提出する。QfO 2022 は 78 種（真核 48、細菌 23、古細菌 7）を含む。共通入力により、アノテーション版、isoform 選択、ID mapping の差を除去する。

### 現行の主な試験

| 試験 | 評価対象 | precision proxy | recall proxy |
|------|----------|-----------------|--------------|
| Species Tree Discordance | sampled ortholog set から再構成した gene tree | species tree との normalized RF distance、incorrect tree fraction | 完成した tree sample 数 / 50,000 試行、予測 pair 数 |
| Generalized Species Tree Discordance（GSTD） | LUCA、Eukaryota、Vertebrata、Fungi 等 | species tree との normalized RF distance | 完成 sample 数または ortholog 数 |
| SwissTree agreement | 人手 curated gene tree が含意する pair | PPV | TPR |
| TreeFam-A agreement | curated/reconciled metazoan gene tree | PPV | TPR |
| VGNC agreement | 人手 curated vertebrate ortholog set | PPV | TPR |
| GO conservation | GO annotation の意味的類似 | 平均 Schlicker similarity | annotation 付き予測数 |
| EC conservation | Enzyme Classification の意味的類似 | 平均 Schlicker similarity | annotation 付き予測数 |
| FAS（2024 追加） | protein feature architecture | 平均 bi-directional FAS | 予測 pair 数 |

QfO のウェブ表示は precision–recall proxy と **Pareto frontier**を用いる。全試験を一つの順位へ潰さず、利用目的ごとに局所的な最適性を比較する。

### 長所

- コミュニティ標準であり、OrthoFinder、SonicParanoid2、FastOMA との直接比較が可能。
- 系統・キュレーション・機能の複数 proxy を同じ入力で評価できる。
- 生データ、既存ツール結果、provenance が OpenEBench にある。

### 限界

- 現行試験の中心は**pairwise relation**であり、orthogroup/HOG の階層構造を十分に評価しない。
- RF distance は、orthology inference の誤りだけでなく、alignment・gene-tree inference error、ILS、HGT にも影響される。
- GO/EC は annotation density と annotation propagation に依存し、orthology prediction を使った注釈による循環性がある。
- curated tree/set は高品質だが、taxon・family coverage が狭い。
- 公開結果は opt-in であり、publication bias がある。
- 2022 set は植物 WGD の代表性が限定的で、OrthoPaiR の主要用途を十分に覆わない。

### OrthoPaiR への採用

**必須。** QfO 2022 Reference Proteomes を固定し、少なくとも全 pair-based tests に提出する。OrthoPaiR が genomic coordinate / GFF を必要とする場合、QfO の locus coordinate と元アノテーションの対応を version-lock して補う。

---

## 3.2 人手 curated ortholog pair（VGNC）

### 方法

VGNC が専門家により curated した vertebrate ortholog set を真値として、予測 pair の PPV と TPR を計算する。2022 QfO 論文では 8 哺乳類、18,282 ortholog sets を収載し、QfO に含まれる 4 種を評価に利用した。olfactory receptor、Cytochrome P450 など複雑な履歴の family も含む。

### 長所

- 自動 orthology tool を根拠にしない人手 curated relation で、循環性が比較的小さい。
- pairwise precision / recall を通常の定義で計算できる。

### 限界

- 主に近縁哺乳類であり、植物・WGD・遠縁比較への外挿はできない。
- 「curated されていない pair」は必ずしも negative ではない。
- 一部のよく研究された family に偏る。

### OrthoPaiR への採用

QfO の一部として採用するが、**植物での優位性の根拠には単独で使わない**。

---

## 3.3 curated reference gene tree（SwissTree / TreeFam-A）

### 方法

専門家が alignment、tree、root、duplication/speciation event を確認した gene tree から ortholog pair を導出し、予測との TP / FP / FN を比較する。

### 指標

- `PPV = TP / (TP + FP)`
- `TPR = TP / (TP + FN)`
- F1、relation class 別（1:1、1:M、M:M）の precision / recall

### 長所

- orthology の進化的定義（speciation node）を直接評価する。
- 特定 family の duplication/loss と difficult cases を解剖できる。

### 限界

- curated tree 自体も tree inference・rooting・taxon sampling の仮説に依存する。
- SwissTree は family 数が少なく、TreeFam-A は metazoan に偏る。
- reference tree を作った方法と評価対象ツールが同系統だと独立性が下がる。

### OrthoPaiR への採用

QfO 標準に加え、植物 family の**独立 gene-tree holdout**を作る。OrthoFinder の tree を OrthoFinder 自身の正解として使うことは避け、可能なら curated trees、synteny、outgroup を併用する。

---

## 3.4 Species Tree Discordance

### 方法

予測 ortholog から 1 種 1 配列をサンプリングし gene tree を再構築する。正しい ortholog set なら、gene tree は既知 species tree と概ね一致するという仮定で RF distance を測る。

### 指標

- normalized Robinson–Foulds distance
- incorrect tree fraction
- completed tree samples / 50,000 trials
- taxonomic level 別の値

### 長所

- 大規模に適用でき、完全な pair-level truth を必要としない。
- precision–coverage trade-off を評価できる。
- QfO の公開値と直接比較可能。

### 限界

- ILS、HGT、hybridization、gene-tree error は「orthology error」でなくても discordance を作る。
- divergent sequence を多く回収する高-recall 法が tree inference error のため不当に罰せられ得る。
- 1 種 1 遺伝子への縮約は co-ortholog/HOG の評価に弱い。

### OrthoPaiR への採用

採用する。ただし植物では accepted species tree が網状進化・WGD の影響を受けるため、RF distance を単独の precision と呼ばず **topological concordance proxy**と記述する。

---

## 3.5 expert-curated orthogroup：Open OrthoBench

### 方法

70 個の Bilateria-level reference orthogroups（RefOGs）と 12 proteomes を使用する。各 RefOG は広い homolog 検索、MSA、rooted gene tree、outgroup を用いて expert curation され、入力 proteomes と評価スクリプトが公開されている。

予測 orthogroup 内の gene pair と RefOG 内 pair を比較する。RefOG サイズによる偏りを避けるため、TP / FP / FN を各 RefOG の `(n - 1)` で正規化する。

### 指標

- normalized pairwise precision、recall、F1
- exact orthogroup recovery 数
- missing genes / missing RefOGs
- **fission error**（一つの RefOG を複数 group に分断）
- **fusion error**（複数 RefOG を一つに混合）
- fragmentation / entropy
- adjusted Rand index
- family size、進化速度、MSA quality、domain complexity 別性能

### 長所

- pair だけでなく group boundary の誤りを測れる。
- 生物学的・技術的 difficulty ごとの stratification が可能。
- FastOMA（2025）、OrthoFinder（2026）でも採用される現行標準。

### 限界

- 70 RefOG、Bilateria 12 種に限定。
- ortholog pair と orthogroup は同じ課題ではない。
- curation の更新で真値が変わる。2020 再解析では旧 RefOG の 44% が修正された。

### OrthoPaiR への採用

`reorgOrthopairs()` の group 出力を評価するため必須。ただし OrthoPaiR は GFF/synteny を使うため、12 species の座標情報を取得できるか先に確認する。

---

## 3.6 機能保存：GO / EC

### 方法

予測 ortholog pair の GO term または EC number の意味的類似度を測り、平均 Schlicker similarity と annotation 付き relation 数を比較する。

### 長所

- 「機能移植に使える ortholog」という利用目的に直結する。
- 多数の family に適用できる。

### 限界

- ortholog でも機能分化し得るため、低い機能類似度は必ずしも誤推定でない。
- paralog でも機能が保存され、高スコアになり得る。
- orthology に基づく electronic annotation propagation が循環性を作る。
- annotation-rich な model organism に偏る。

### OrthoPaiR への採用

補助評価として採用する。experimental evidence code のみの解析と全 GO の解析を分け、注釈伝播による循環を明示する。

---

## 3.7 Feature Architecture Similarity（FAS）

### 方法

Pfam / SMART domain、signal peptide、transmembrane region、low-complexity region 等を多層 feature architecture として表し、重複 annotation を DAG 上の score maximization で解決する。方向付き score を双方向平均し、0（共有 feature なし）から 1（同一または包含 architecture）で評価する。

### 指標

- pair ごとの bidirectional FAS
- ツールごとの mean / median FAS
- predicted pair 数との Pareto 関係
- relation class・系統距離・protein length ratio 別 FAS

### 長所

- GO/EC より annotation propagation の循環が小さい。
- domain gain/loss、fusion/fission、partial orthology を検出しやすい。
- QfO 2024 の標準試験で競合ツール値がある。

### 限界

- 高 FAS は orthology を証明しない。paralog でも architecture は保存される。
- feature database の coverage と gene model quality に依存する。
- 真の domain innovation を誤りとして罰し得る。

### OrthoPaiR への採用

必須の補助評価。特に `Only_OrthoFinder` / `Only_Sonic` / `Only_OrthoPaiR` の差分集合で FAS 分布を比較する。

---

## 3.8 シミュレーションによる既知真値

### 方法

ALF 等で既知 species tree 上に配列を進化させ、duplication、loss、HGT、substitution rate、rate heterogeneity、indel、domain event、WGD、gene model error を制御する。真の gene tree と speciation/duplication event が既知なので、全 relation の TP / FP / FN を測れる。

Natsidis et al.（2021）は 57 metazoans、5,000 genes、200 simulations を用い、loss/duplication を起こさない単純設定でも、進化速度と site-rate heterogeneity により systematic orthology error が生じることを示した。

### 指標

- pairwise precision、recall、F1、MCC
- 1:1 / 1:M / M:M relation class 別性能
- orthogroup exact recovery、fission、fusion
- duplication / loss event の precision / recall
- branch length、evolutionary rate、WGD 後時間、欠損率別性能
- downstream gain/loss error、phylostratigraphic age error

### 長所

- 唯一、全 relation の完全な真値を持てる。
- 失敗条件を一因子ずつ変えられる。
- algorithm tuning と regression test に向く。

### 限界

- simulator と parameter が現実を再現しない場合、simulation に過適合する。
- synteny、chromosome rearrangement、annotation error、isoform complexity を十分再現しないモデルが多い。
- simulator の仮定と同じ仮定を持つツールが有利になる。

### OrthoPaiR への採用

必須。次の直交する系列を作る。

1. **単一コピー・速度勾配**
2. **small-scale duplication/loss**
3. **WGD + differential loss（pseudoortholog を明示的に生成）**
4. **染色体 rearrangement と synteny decay**
5. **gene model fragmentation / isoform 違い**
6. **taxon sampling の疎密**

シミュレーション生成コードと seed を固定し、開発 set と最終 holdout set を分ける。

### 補足：reconciliation / ILS シミュレーション

QfO の STD は「正しい ortholog set なら gene tree ≈ species tree」を仮定する。ILS、HGT、gene-tree error がある場合、この仮定は破れる。そのため 2021–2026 年には、**種系統再構成・event recovery** を主目的としたシミュレーション評価も並行して発展した。

| 枠組み | 文献 | 測るもの | OrthoPaiR への位置づけ |
|--------|------|----------|------------------------|
| ALF | Natsidis et al. 2021 | 真の pair/group に対する PR；下流 gain/loss への影響 | pair/group 精度の主シミュレーション |
| Zombi | 近年の phylogenomics 論文で利用 | DTL 付き genome/family evolution | WGD 以外の複雑シナリオの補助 |
| SpeciesRax / AleRax | Morel et al. 2022 / 2024 | DTL 下での種系統 RF、reconciliation 精度 | reorg の種系統一貫性オプション評価に有用 |
| reconcILS / dupcoal | Mishra et al. 2023– | ILS と DL event の回収 | STD RF を過大解釈しないための cautionary benchmark |

**要点:** simulation はコミュニティ標準の bake-off としては QfO ほど整備されていないが、**mechanistic error study** と **pseudoortholog / ILS の分離**には不可欠である。

---

## 3.9 シンテニー・倍数性・遺伝子順序による検証

### 方法

近縁ゲノム、倍数体、既知 WGD 系統で以下を測る。

- conserved collinear block 内の ortholog 率
- predicted pair の近傍保存率
- single-copy group から作った gene tree の species topology 一致率
- allopolyploid を一種として扱う場合と、subgenome に分割して別種扱いする場合の group 一致率
- expected copy-number / ploidy ratio からの逸脱

pSONIC（2021）は cotton tribe で、single-copy group の約 99% が期待 species topology に一致し、allopolyploid の subgenome 分割前後で 97.5% の group が一致したと報告した。

### 長所

- 植物 WGD と pseudoortholog に直接対応する。
- 配列類似度・gene tree とは独立性のある genomic-context 証拠。
- OrthoPaiR の設計上の強みを検証できる。

### 限界

- synteny は orthology の必要条件ではなく、遠縁・高 rearrangement 系統では崩れる。
- tandem duplication、assembly error、annotation error に影響される。
- OrthoPaiR が synteny を使用するため、同じ synteny で評価すると循環性が生じる。

### OrthoPaiR への採用

必須だが、**推定に使った anchor/block と同じデータだけで評価しない**。独立手法の collinearity、curated WGD history、gene tree、outgroup を組み合わせる。

---

## 3.10 ツール間一致・不一致ネットワーク

### 方法

全ツールの pair 集合から consensus、tool-specific、disputed relations を作る。

```text
Core       = 全ツール共通
Only_A     = A \ (他ツールの和集合)
Disputed   = 対象ツール間の対称差
Support_k  = k 個以上のツールが支持
```

QfO 2024 は 18 手法の pair overlap、unique pair 数、pair-vector / benchmark-vector の PCA を評価した。Yang et al.（2025）は genome context の unsupervised clustering から Signal Jaccard Index（SJI）を構成し、protein network の degree centrality が ortholog database 間の不一致と関連することを示した。

### 指標

- Jaccard / asymmetric overlap
- support count
- unique pair fraction
- pairwise agreement matrix
- SJI、degree centrality
- disputed pair の curated-tree / FAS / synteny 支持率

### 長所

- ground truth がなくても、どの relation が不安定かを発見できる。
- method-specific bias と family-specific difficulty を可視化できる。

### 限界

- 多数決は真値ではない。類似アルゴリズム群が同じ誤りを共有する。
- 同じ BLAST/DIAMOND、MCL、database annotation を使うツール間は独立でない。
- consensus を正解にすると、新規性のある正しい予測を罰する。

### OrthoPaiR への採用

探索・エラー解剖には採用するが、主要 accuracy claim には単独使用しない。`Only_OF`、`Only_Sonic`、`Only_OP` を Tier B/C の独立証拠で評価する。

---

## 3.11 下流解析への影響

### 方法

ortholog/orthogroup の出力を実際の下流解析へ投入し、結果が既知の biology または simulation truth を再現するか測る。

- species tree topology / branch support
- gene gain/loss 数
- phylostratigraphic gene age
- ancestral gene content
- functional annotation transfer accuracy
- conserved single-copy marker の回収

Natsidis et al.（2021）は、真の gain/loss が無い simulation でも orthology inference error だけで見かけの gain/loss と phylostratigraphic pattern が生じ、経験データに似た signal を作り得ることを示した。

### 長所

- ユーザーが最終的に必要とする解析品質を評価できる。
- 小さな pair-level error の実質的影響を測れる。

### 限界

- 下流モデル自身の誤差と orthology error を分離しにくい。
- 一つの用途に良い tool が別用途にも良いとは限らない。

### OrthoPaiR への採用

植物の comparative genomics 用途に合わせ、少なくとも species tree、gain/loss、single-copy marker の三つを評価する。

---

## 3.12 計算性能・スケーラビリティ

### 方法

同一 hardware、同一 container、同一 input、同一 thread 数で wall time と resource を測る。SonicParanoid2（2024）は QfO、200 eukaryotes、2,000 MAGs を用い、FastOMA（2025）は 2,086 eukaryotic proteomes まで評価した。

### 必須指標

| 指標 | 備考 |
|------|------|
| wall-clock time | end-to-end と各 phase の両方 |
| CPU time / core-hours | thread 数が違う比較に必要 |
| peak RSS | child process を含める |
| disk peak / I/O bytes | all-vs-all 中間ファイルを含める |
| completion rate | timeout / OOM を失敗として含める |
| scaling exponent | 種数・protein 数を増やした log–log slope |
| output size | pair 数、group 数、relation class |
| accuracy per cost | F1/core-hour、TP/GiB 等。ただし単一総合点だけにしない |

### 注意

- similarity search の既計算結果を再利用する場合と end-to-end を分ける。
- database cache の warm/cold を区別する。
- default、fast、sensitive mode は別条件として扱う。
- tool ごとに異なる出力（pair、OG、HOG）を同じ「件数」で比較しない。

---

## 4. ベンチマーク設計上の共通注意

### 4.1 出力単位を混ぜない

| 出力 | 評価単位 |
|------|----------|
| pairwise ortholog | unordered cross-species gene pair |
| co-ortholog | 1:M / M:M relation を明示 |
| flat orthogroup | group membership、fission/fusion |
| HOG | taxonomic level ごとの nested group |
| gene tree | duplication/speciation event、tree topology |

OG を全直積 pair に展開すると、大きな group が FP と pair 数を支配する。pairwise benchmark と group benchmark を別々に報告する。

### 4.2 「precision proxy」を precision と断定しない

Species-tree RF、GO/EC、FAS、synteny は真の pair-level precision ではない。報告ではそれぞれ、

- topological concordance
- functional similarity
- feature-architecture similarity
- genomic-context support

と明記する。

### 4.3 開発 set と holdout set を分離する

同じ BUSCO/VGNC/plant pairs を見ながら閾値を調整すると benchmark overfitting になる。

1. **training/development:** parameter sweep と FN audit に利用
2. **validation:** algorithm family の選択
3. **locked holdout:** 最終比較で一度だけ利用

QfO public tests も繰り返し最適化すると事実上 development set になるため、独立 simulation seeds と独立 plant clades を holdout にする。

### 4.4 共通入力と provenance を固定する

- genome/proteome assembly accession
- annotation release
- canonical isoform rule
- gene ↔ transcript ↔ protein mapping
- tool version / commit hash
- container digest
- parameters、seed、thread 数
- species tree
- hardware / OS
- timeout / memory limit

### 4.5 信頼区間と統計単位

pair を独立標本とみなすと、巨大 gene family が過大な重みを持つ。gene family または species pair を cluster とした bootstrap を使い、95% CI を報告する。多重比較時は FDR を調整する。

---

## 5. OrthoPaiR 改修前に推奨するベンチマーク構成

## 5.1 最小構成

| 層 | データ | 主指標 | 目的 |
|----|--------|--------|------|
| A | QfO 2022 | 全標準試験の Pareto 位置 | 国際標準・競合との直接比較 |
| B | Open OrthoBench | precision/recall/F1、exact、fission/fusion | reorg group boundary |
| C | 植物 3–4 clade | curated/synteny/tree-supported pair | WGD・遠縁・pseudoortholog |
| D | ALF/独自 simulation | 真の pair/group/event 精度 | 完全真値と failure mode |
| E | resource scaling | time/RSS/disk/completion | 実用性 |

## 5.2 植物特化セット

最低でも以下を別々に用意する。

1. **近縁 diploid pair** — synteny が強く、1:1 が多い
2. **WGD 後の diploidized pair** — differential loss / pseudoortholog
3. **allopolyploid + progenitors + outgroup** — ploidy-aware consistency
4. **遠縁 pair** — synteny decay と sequence divergence
5. **annotation-quality perturbation** — gene fragmentation、missing gene、代表 transcript 変更

### 正例

- curated gene tree の speciation relation
- high-confidence synteny + outgroup-supported relation
- QfO/VGNC 相当の専門家 curated plant family
- simulation truth

### 負例

- duplication node を挟む paralog/pseudoortholog
- WGD の別 subgenome copy
- gene tree + synteny + outgroup が一致する curated negative

「正例リストに無いもの」をすべて負例としない。negative は明示的に curated する。

## 5.3 比較対象

- OrthoPaiR pairwise
- OrthoPaiR reorg
- OrthoFinder（固定版）
- SonicParanoid2（fast/default/sensitive を区別）
- FastOMA
- OMA Groups / HOGs
- BBH/RBH baseline

必要に応じて Proteinortho、Broccoli、OrthoMCL を追加する。

## 5.4 必須 ablation

| Ablation | 確認する問い |
|----------|--------------|
| pairwise vs reorg | 多種グラフは本当に改善するか |
| synteny/anchor on vs off | genomic context の純効果 |
| bridge inference on vs off | recall gain と FP 増加 |
| graph_qc on vs off | pseudoortholog 抑制効果 |
| strict / broad / all | confidence label が実精度を層別化できるか |
| protein-only vs GFF-assisted | 入力情報量の差を明示 |

## 5.5 主要評価表

各 dataset × relation class について以下を出す。

- TP / FP / FN / TN（TN が定義可能な場合のみ）
- precision / recall / F1 / MCC
- PR curve と Pareto frontier
- pair 数 / orthogroup 数
- exact group、fission、fusion、entropy
- tree concordance、FAS、synteny support
- wall time、core-hours、peak RSS、disk
- family-level bootstrap 95% CI

---

## 6. 実施順序

### Phase 0 — 評価仕様の固定

- [ ] pair / OG / HOG の評価単位を別々に定義
- [ ] primary endpoint と secondary endpoint を事前登録
- [ ] development / validation / holdout を分割
- [ ] competitor version、container、hardware を固定
- [ ] pair ID 正規化と isoform policy を固定

### Phase 1 — 外部標準

- [ ] QfO 2022 Reference Proteomes 上で OrthoPaiR を実行
- [ ] OrthoXML または pair TSV を OpenEBench に提出
- [ ] 8 種類の標準 benchmark 結果を保存
- [ ] Open OrthoBench を reorg group 出力で実行

### Phase 2 — 植物 benchmark

- [ ] WGD/pseudoortholog を含む curated positive/negative set を作成
- [ ] synteny と gene tree を独立に作成
- [ ] 近縁・遠縁・倍数体を層別化
- [ ] OrthoFinder / SonicParanoid2 / FastOMA と共通入力で比較

### Phase 3 — simulation と robustness

- [ ] ALF/代替 simulator の parameter grid
- [ ] WGD + differential loss simulation
- [ ] rearrangement / annotation error perturbation
- [ ] species sampling と genome 数の sensitivity

### Phase 4 — 最終 holdout

- [ ] algorithm と threshold を freeze
- [ ] holdout を一度だけ評価
- [ ] CI、effect size、failure cases を報告
- [ ] 全入力・予測・評価コードを archive

---

## 7. OrthoPaiR で最優先すべき判断

アルゴリズム改修に戻る前に、次を決める必要がある。

1. **主対象は pairwise ortholog か、orthogroup/HOG か。**
2. **狭義 1:1 precision と co-ortholog recall のどちらを primary endpoint にするか。**
3. **植物 WGD/pseudoortholog をどの curated set または simulation で評価するか。**
4. **QfO を開発指標にするか、locked external test にするか。**
5. **競合への優位性を accuracy、pair-efficiency、runtime のどの Pareto 面で主張するか。**

推奨 primary claim は次の形である。

> OrthoPaiR は、植物の WGD・遠縁比較において、同等の recall 水準で pseudoortholog precision、group fusion 抑制、synteny/tree support、および計算資源の Pareto frontier に位置する。

この主張を支持できる評価系を先に確立し、その後に link_all gap や reorg heuristic を改修する。

---

## 8. 参考 URL

- [QfO Orthology Benchmarking Service](https://orthology.benchmarkservice.org/)
- [QfO documentation](https://orthology.benchmarkservice.org/proxy/doc)
- [QfO 2022 public project](https://orthology.benchmarkservice.org/proxy/projects/2022/)
- [Open OrthoBench](https://github.com/davidemms/Open_Orthobench)
- [FAS](https://github.com/BIONF/FAS)
- [pSONIC](https://github.com/conJUSTover/pSONIC)

---

## 9. 調査上の留保

- 文献検索日は 2026-07-18。2026 年論文は early online を含む。
- QfO のウェブ UI では taxonomic subset ごとに複数試験が展開されるため、「7/8 種類」「11/12 tests」など論文間で数え方が異なる。本書では**方法の種類**として、従来 7 種類に 2024 年追加の FAS を加えた 8 種類で整理した。
- QfO の各指標は universal ground truth ではなく proxy である。単一の総合順位を作らず、precision–recall–cost の Pareto 比較を基本とする。
