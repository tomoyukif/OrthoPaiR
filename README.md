## OrthoPaiR

[![R](https://img.shields.io/badge/R-%3E%3D4.1-blue)](https://www.r-project.org/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

**OrthoPaiR** is an R package for **synteny‑aware ortholog detection** across two or more genomes.
The current workflow is file-based under a `working_dir` and consists of:
input reformatting, reciprocal search (RBH), pairwise ortholog inference, and optional multi-genome reorganisation.

### Features

- **Multi‑genome support**: run ortholog detection for all pairwise combinations of genomes.
- **Synteny‑aware pairing**: combine reciprocal best hits with genomic context.
- **Pairwise ortholog tables**: writes per-pair TSV files in `working_dir/orthopair`.
- **Network‑based ortholog groups**: reorganisation step writes graph and orthogroup tables in `working_dir/reorg_out`.

### Prerequisites

- **R packages** (installed via CRAN and Bioconductor):

```r
install.packages("devtools")
install.packages("BiocManager")

BiocManager::install(c(
  "GenomicRanges", "GenomeInfoDb", "GenomicFeatures",
  "BiocGenerics", "Biobase", "IRanges",
  "rtracklayer", "Biostrings", "S4Vectors", "BSgenome"
))

install.packages(c("parallel", "ggplot2", "dplyr"))
```

- **External tools**
  - [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs) (`makeblastdb`, `blastn` in `blast_path`) – required for DNA sequence similarity; `blastp` is required only when `use_prot = TRUE` for protein sequence similarity
  - [DIAMOND](https://github.com/bbuchfink/diamond) (`diamond` in `diamond_path`) – required when `use_prot = TRUE` for protein sequence similarity
  - [Miniprot](https://github.com/lh3/miniprot) (optional; only if `run_miniprot = TRUE` to predict missing ORFs) – requires genome FASTA files in `orgInputFiles()`

### Installation

Install the development version from GitHub:

```r
install.packages("devtools")
devtools::install_github("tomoyukif/OrthoPaiR", dependencies = TRUE)
```

### Usage (simplified)

The main workflow is:

1. Prepare an input data.frame (`name`, `genome`, `gff`, `cds`, `prot`).
2. Reformat and index inputs with `reformatFiles()`.
3. Run reciprocal search with `rbh()`.
4. Infer pairwise orthologs with `orthopair()`.
5. Optionally reorganise all-pair outputs with `reorgOrthopairs()`.

```r
library(OrthoPaiR)

wd <- "path/to/output/orthopair"
input_df <- data.frame(
  name = c("Osat", "Hvul"),
  genome = c("path/to/osat_genome.fa", "path/to/hvul_genome.fa"),
  gff = c("path/to/osat.gff3", "path/to/hvul.gff3"),
  cds = c("path/to/osat_cds.fa", "path/to/hvul_cds.fa"),
  prot = c(NA, NA),
  stringsAsFactors = FALSE
)

reformatFiles(
  object = input_df,
  working_dir = wd,
  overwrite = TRUE,
  n_threads = 8
)

rbh(
  working_dir = wd,
  blast_path = "/path/to/blast/bin",
  target_pair = NULL,
  n_threads = 8,
  overwrite = TRUE
)

orthopair(
  working_dir = wd,
  n_threads = 8
)

reorgOrthopairs(
  working_dir = wd,
  rename = TRUE,
  n_threads = 8
)
```

- For **two genomes**, `orthopair()` output in `working_dir/orthopair/<genomeA>_<genomeB>.tsv` is usually sufficient.
- For **three or more genomes**, run `reorgOrthopairs()` to create graph and orthogroup summaries in `working_dir/reorg_out`.

For a complete four-genome example, refer to the main vignette:

```r
vignette("orthopair", package = "OrthoPaiR")
```

When you run the workflow, the output directory will have the following structure:

```
/your_working_dir/
├── blast/
│   ├── *_group_cds.fa
│   ├── *.blastdb*
│   ├── *_blast.out
│   └── ...
├── rbh/
│   ├── all.rbh.tsv
│   ├── 1001_1002.rbh
│   └── ...
├── orthopair/
│   ├── 1001_1002.tsv
│   ├── 1001_1003.tsv
│   ├── orthopair_pairwise_mutual_ci_stats.tsv
│   └── orthopair_genome_mean_mutual_ci_matrix.tsv
├── reorg_out/
│   ├── orthopair.graphml
│   ├── orthopair_list.tsv
│   └── pairwise/
│       ├── 1001_1002.tsv
│       └── ...
├── input/
│   ├── 1_Osat/
│   ├── 2_Hvul/
│   └── ...
```

### Documentation

For a detailed description of input requirements, pipeline steps, and all output files, see the vignette:

```r
vignette("orthopair", package = "OrthoPaiR")
```

### License

OrthoPaiR is licensed under the **GPL‑3.0** License.  
See `LICENSE.md` for details.

